function [Accuracy, Avg_F_measure, Norm_F] = Results_statistics (label, Predicted)

error = 0;
for i=1: length(label)
    if label(i) ~= Predicted(i)
        error = error+1;    
    end;
end;
error_rate = error/length(label);
Accuracy   = 1 - error_rate;

[uniques,numUnique] = count_unique(label);
[uniques_pre,numUnique_Pre] = count_unique(Predicted);

if length(uniques) > length(uniques_pre)
    disp('Predicted has less number of Classes ');
end;
if length(uniques) < length(uniques_pre)
    disp('Predicted has more number of Classes ');
end;
Avg_accuracy  = 0;
Avg_F_measure = 0;
Accuracies    = zeros(length(uniques),1);
F_measures    = zeros(length(uniques),1);
for i =1: length(uniques)
    
    label_index = find(label==uniques(i)); %actual label
    Pre_index   = find(Predicted==uniques(i));  %predicted label
    
    True_possitive = length( intersect(label_index,Pre_index) );
    False_Negative = length(label_index)-True_possitive;
    False_Positive = length( setdiff(Pre_index, label_index) );
    
    Precision      = True_possitive/length(Pre_index);
    Recall         = True_possitive/length(label_index);
    
    F_Measure      = 2*Precision*Recall/(Precision+Recall);
    if isnan(F_Measure)
        F_Measure = 0;
    end;
    
    accuracy       = True_possitive/length(label_index);
    Avg_accuracy = Avg_accuracy +  accuracy;
    Avg_F_measure = Avg_F_measure + F_Measure;
    
    Accuracies(i) = accuracy;
    F_measures(i) = F_Measure;
    
    class_errors{i}.label = uniques(i);
    class_errors{i}.instance_no = length(label_index);
    class_errors{i}.True_possitive = True_possitive;
    class_errors{i}.False_Negative = False_Negative;
    class_errors{i}.False_Positive = False_Positive;
    class_errors{i}.Precision = Precision;
    class_errors{i}.Recall = Recall;
    class_errors{i}.accuracy = accuracy;
    class_errors{i}.F_Measure = F_Measure;
end;

Avg_accuracy = Avg_accuracy/length(uniques);
Avg_F_measure = Avg_F_measure /length(uniques);

a = F_measures.*numUnique;
Norm_F = sum(a)/sum(numUnique);

%% compute confusion matrix
[label_uni,label_inst] = count_unique(label);
[pre_uni,pre_inst] = count_unique(Predicted);
Confus_matrix = zeros(length(label_uni),length(pre_uni));
for i=1: length(label_uni)
    for j=1: length(pre_uni)
        label_Index = find(label ==label_uni(i));
        pre_Index = find(Predicted ==label_uni(j));
        temp = length( intersect(label_Index,pre_Index));
        Confus_matrix(i,j)= temp/label_inst(i)*100;
        
    end;
end;    
    
    

end


function [uniques,numUnique] = count_unique(x,option)
%COUNT_UNIQUE  Determines unique values, and counts occurrences
%   [uniques,numUnique] = count_unique(x)
%
%   This function determines unique values of an array, and also counts the
%   number of instances of those values.
%
%   This uses the MATLAB builtin function accumarray, and is faster than
%   MATLAB's unique function for intermediate to large sizes of arrays for integer values.  
%   Unlike 'unique' it cannot be used to determine if rows are unique or 
%   operate on cell arrays.
%
%   If float values are passed, it uses MATLAB's logic builtin unique function to
%   determine unique values, and then to count instances.
%
%   Descriptions of Input Variables:
%   x:  Input vector or matrix, N-D.  Must be a type acceptable to
%       accumarray, numeric, logical, char, scalar, or cell array of
%       strings.
%   option: Acceptable values currently only 'float'.  If 'float' is
%           specified, the input x vector will be treated as containing
%           decimal values, regardless of whether it is a float array type.
%
%   Descriptions of Output Variables:
%   uniques:    sorted unique values
%   numUnique:  number of instances of each unique value
%
%   Example(s):
%   >> [uniques] = count_unique(largeArray);
%   >> [uniques,numUnique] = count_unique(largeArray);
%
%   See also: unique, accumarray

% Author: Anthony Kendall
% Contact: anthony [dot] kendall [at] gmail [dot] com
% Created: 2009-03-17

testFloat = false;
if nargin == 2 && strcmpi(option,'float')
    testFloat = true;
end

nOut = nargout;
if testFloat
    if nOut < 2
        [uniques] = float_cell_unique(x,nOut);
    else
        [uniques,numUnique] = float_cell_unique(x,nOut);
    end
else
    try %this will fail if the array is float or cell
        if nOut < 2
            [uniques] = int_log_unique(x,nOut);
        else
            [uniques,numUnique] = int_log_unique(x,nOut);
        end
    catch %default to standard approach
        if nOut < 2
            [uniques] = float_cell_unique(x,nOut);
        else
            [uniques,numUnique] = float_cell_unique(x,nOut);
        end
    end
end

end

function [uniques,numUnique] = int_log_unique(x,nOut)
%First, determine the offset for negative values
minVal = min(x(:));

%Check to see if accumarray is appropriate for this function
maxIndex = max(x(:)) - minVal + 1;
if maxIndex / numel(x) > 1000
    error('Accumarray is inefficient for arrays when ind values are >> than the number of elements')
end

%Now, offset to get the index
index = x(:) - minVal + 1;

%Count the occurrences of each index value
numUnique = accumarray(index,1);

%Get the values which occur more than once
uniqueInd = (1:length(numUnique))';
uniques = uniqueInd(numUnique>0) + minVal - 1;

if nOut == 2
    %Trim the numUnique array
    numUnique = numUnique(numUnique>0);
end
end 

function [uniques,numUnique] = float_cell_unique(x,nOut)

if ~iscell(x)
    %First, sort the input vector
    x = sort(x(:));
    numelX = numel(x);
    
    %Check to see if the array type needs to be converted to double
    currClass = class(x);
    isdouble = strcmp(currClass,'double');
    
    if ~isdouble
        x = double(x);
    end
    
    %Check to see if there are any NaNs or Infs, sort returns these either at
    %the beginning or end of an array
    if isnan(x(1)) || isinf(x(1)) || isnan(x(numelX)) || isinf(x(numelX))
        %Check to see if the array contains nans or infs
        xnan = isnan(x);
        xinf = isinf(x);
        testRep = xnan | xinf;
        
        %Remove all of these from the array
        x = x(~testRep);
    end
    
    %Determine break locations of unique values
    uniqueLocs = [true;diff(x) ~= 0];
else
    isdouble = true; %just to avoid conversion on finish
    
    %Sort the rows of the cell array
    x = sort(x(:));
    
    %Determine unique location values
    uniqueLocs = [true;~strcmp(x(1:end-1),x(2:end)) ~= 0] ;
end

%Determine the unique values
uniques = x(uniqueLocs);

if ~isdouble
    x = feval(currClass,x);
end

%Count the number of duplicate values
if nOut == 2
    numUnique = diff([find(uniqueLocs);length(x)+1]);
end
end
