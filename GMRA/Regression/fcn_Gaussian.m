function Y = fcn_Gaussian( X )

R = 0.1;

idx = randperm(size(X,2),1);
distsq = sum(bsxfun(@minus,X,X(:,idx)).^2,1);

Y = exp(-distsq/(R*EstimateSetDiam(X))^2)';

return