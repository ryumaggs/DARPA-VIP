function [medianmindist,medianmediandist,minmindist] = EstimateMinPairwiseDistance( X )

%
% function medianmindist = EstimateMinPairwiseDistance( X )
%
% Simple cheap randomized estimator for the median of the minimum distance between a point and its nearest neighbor
%

NPTSTOSAMPLE        = 100;

idxs = randperm( size(X,2), min(NPTSTOSAMPLE,size(X,2)) );

Xsub = X(:,idxs);

diamsq = 0;

for k = 1:NPTSTOSAMPLE
    distsq = sum(bsxfun(@minus,X,Xsub(:,k)).^2,1);
    mediandistsq(k) = median(distsq);
    mindistsq(k)    = min(distsq(distsq>0));    
end

medianmindist       = median(sqrt(mindistsq));
medianmediandist    = median(sqrt(mediandistsq));
minmindist          = sqrt(min(mindistsq));

return