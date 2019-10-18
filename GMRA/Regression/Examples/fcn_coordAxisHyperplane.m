function Y = fcn_coordAxisHyperplane( X )

MAX_TRIALS = 100;

ntrials = 0;
while true
    delta = zeros(1,size(X,1));
    delta(randperm(size(X,1),1))=1;
    
    Y = delta*bsxfun(@minus,X,mean(X,2))/EstimateSetDiam(X);
    if max(Y)-min(Y)>10*eps
        break;
    end
    ntrials = ntrials + 1;
    if ntrials > MAX_TRIALS
        warning('fcn_coordAxisHyperplane: constant function');
        break;
    end
end

Y = Y';

return