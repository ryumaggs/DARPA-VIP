function Y = fcn_oscpole( X )

DISTANCEFROMMANIFOLD = 0.02;
OSCILLATIONRATE      = 5;

idx = randperm( size(X,2),1 );

delta = randn(size(X,1),1);
delta = delta/norm(delta);

pole = X(:,idx) + DISTANCEFROMMANIFOLD*delta;

Y = sqrt(sum(bsxfun(@minus,X,pole).^2,1))';

Y = 1./(Y.^2).*sin(2*pi*OSCILLATIONRATE*Y/(0.5*EstimateSetDiam(X)));

return