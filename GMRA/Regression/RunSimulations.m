

N = 10;

for i = 1:N
    RunExamples;
    GMRAu = min(squeeze(unifError.rel(:,:,k)),[],2)';
    GMRAa = min(squeeze(adaptError.rel(:,:,k)),[],2)';
    nn = min(squeeze(Error_NN(k,:)));
    nystrom = min(squeeze(Error_nystromCoRe(k,:)));
    cart = Error_regTree(k);
end