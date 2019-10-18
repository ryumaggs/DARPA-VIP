function Y = fcn_pole( X )

DISTANCEFROMMANIFOLD = 0.1;

idx = randperm( size(X,2),1 );

delta = randn(size(X,1),1);
delta = delta/norm(delta);

pole = X(:,idx) + DISTANCEFROMMANIFOLD*delta;

Y = sqrt(sum(bsxfun(@minus,X,pole).^2,1))';

Y = 1./Y;

return