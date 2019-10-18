function diam = EstimateSetDiam( X )

%
% function diam = EstimateSetDiam( X )
%
% Simple cheap randomized estimator for the diameter of a set
%

NPTSTOSAMPLE        = 1000;
NITERS              = 10;
STOPPINGCRITERION   = 1e-3;

idxs = randperm( size(X,2), min(NPTSTOSAMPLE,size(X,2)) );

Xsub = X(:,idxs);

diamsq = 0;
curidx = 1;

for k = 1:NITERS
    distsq = sum(bsxfun(@minus,Xsub,Xsub(:,curidx)).^2,1);
    [newdiamsq,curidx] = max(distsq);  
    if k>1 && newdiamsq/diamsq<STOPPINGCRITERION
        diamsq = newdiamsq;
        break;
    end
    diamsq = newdiamsq;
end

diam = sqrt(diamsq);

return