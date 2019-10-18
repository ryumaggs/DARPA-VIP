function Error = ComputeRegressionError(Y_test,Yhat,uniqueLabels)

%
% function Error = ComputeRegressionError(Y_test,Yhat,uniqueLabels)
%
% IN:
%   Y_test  : test values
%   Yhat    : estimated test values
%   [uniqueLabels] : should be provided, in case of classification, in order to map Yhat to labels
%
% OUT:
%   Error.abs   : absolute mean L^2 error, #degrees x #partitions x #functions
%   Error.rel   : normalized mean L^ 2 error (relative to L^2 norm of input)

if nargin<3 || isempty(uniqueLabels)
    Error.abs = zeros(length(Yhat.Yhat),max(cellfun(@length,Yhat.Yhat)),size(Y_test,2));
    Error.rel = zeros(length(Yhat.Yhat),max(cellfun(@length,Yhat.Yhat)),size(Y_test,2));
    lYnorms = sqrt(sum(Y_test.^2,1));
    for i = 1:length(Yhat.Yhat)
        for l = 1:length(Yhat.Yhat{i})
            Error.abs(i,l,:) = sqrt(sum((Yhat.Yhat{i}{l}-Y_test).^2,1));
        end
        Error.rel(i,:,:) = bsxfun(@rdivide,squeeze(Error.abs(i,:,:)),lYnorms);
    end
        
    Error.abs = Error.abs/sqrt(size(Y_test,1));
else
    for i = 1:length(Yhat.Yhat)
        
        for l = 1:length(Yhat.Yhat{i})
            Error.abs(i,l) = sum(mapToLabels(Yhat.Yhat{i}{l},uniqueLabels)~=Y_test);
        end        
    end
    Error.abs = Error.abs/length(Y_test);
    Error.rel = Error.abs;    
end

return