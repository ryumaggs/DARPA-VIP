function Y = fcn_norm( X )

Y = sqrt(sum(X.^2,1))';

return