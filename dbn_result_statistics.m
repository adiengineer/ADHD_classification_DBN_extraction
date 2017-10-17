function[Accuracy, Ypre] = dbn_result_statistics(dbn_sizes,train_data, test_data, train_label, test_label, opts, n_out)

%%%Set up, Train DBN and use its weights for top layer NN initialization 
%%%the dbn.sizes is the number of hidden neurons for each RBM layers
%%%and adjustable
dbn.sizes=dbn_sizes
dbn = dbnsetup(dbn, train_data, opts);
dbn = dbntrain(dbn, train_data, opts);

nn = dbnunfoldtonn(dbn, n_out);
nn.activation_function = 'sigm';

nn = nntrain(nn, train_data, train_label, opts);
% [er, bad] = nntest(nn, test_data, test_label);

Ytst_pre = nnpredict(nn, test_data);

% Double the result 
[~, expected] = max(test_label,[],2);
[acc, af, nf] = Results_statistics (expected, Ytst_pre);

Accuracy=acc
Ypre=Ytst_pre

end
