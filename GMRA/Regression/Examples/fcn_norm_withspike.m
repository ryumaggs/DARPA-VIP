function Y = fcn_norm_withspike( X )

pole = [0.9;0;0.5;zeros(size(X,1)-3,1)];

Y = sqrt(sum(bsxfun(@minus,X,pole).^2,1))';

Y = 1./Y;

% Y = sqrt(sum(X.^2,1))';
% Y = 1./(0.1-Y);

return