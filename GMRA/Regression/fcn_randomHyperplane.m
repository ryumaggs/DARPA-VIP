function Y = fcn_randomHyperplane( X )

delta = randn(1,size(X,1));
delta = delta/norm(delta);

Y = delta*bsxfun(@minus,X,mean(X,2))/EstimateSetDiam(X);
Y = Y';

return