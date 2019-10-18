function [Y,meanX,normsX] = Whiten( X )

meanX               = mean(X,1);
Y                   = bsxfun(@minus,X,meanX);
normsX              = sqrt(sum(Y.^2,1));
Y                   = bsxfun(@rdivide, Y, normsX);
Y(:,normsX<eps)     = 0;


return