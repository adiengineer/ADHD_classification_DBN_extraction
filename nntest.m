function [er, bad] = nntest(nn, x, y)
    labels = nnpredict(nn, x);
%     predictions = nnpredict(nn, x);
%     labels = getypre(predictions,n,ny,ylag);
    label_size=size(labels)    
    [~, expected] = max(y,[],2);
    expected_size=size(expected)
    bad = find(labels ~= expected);    
    er = numel(bad) / size(x, 1);
end
