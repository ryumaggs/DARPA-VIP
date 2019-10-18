function Y = fcn_BallIndicator( X )

diamX   = EstimateSetDiam(X);
R       = 0.2;

for loopIdx=1:1000
    idx     = randperm(size(X,2),1);
    distsq  = sum(bsxfun(@minus,X,X(:,idx)).^2,1);
    
    idxs    = find(distsq<=(R*diamX)^2);
    if length(idxs)>0.10*size(X,2),
        break
    end
end

Y       = zeros(size(X,2),1);
Y(idxs) = 1;

return