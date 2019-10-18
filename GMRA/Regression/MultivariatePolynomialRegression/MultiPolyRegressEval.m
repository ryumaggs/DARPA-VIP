function [yhat,NewScores] = MultiPolyRegressEval( poly, X, NewScores )

if isempty( poly.Coefficients ), yhat = []; return; end;
if isempty( X ), yhat = zeros(size(X,2),1); return; end;
if nargin<3, NewScores = []; end

yhat = zeros(size(X,2),size(poly.Coefficients,2));

if poly.Lim>1
    if isempty(NewScores)
        for k = 1:size( X,2 )
            NewScores   = repmat(X(:,k)',[length(poly.PowerMatrix) 1]).^poly.PowerMatrix;
            NewScores   = prod(NewScores,2);
            yhat(k,:)   = NewScores'*poly.Coefficients;                                                                         % The estimate for the new data point
        end
    else
        yhat   = NewScores'*poly.Coefficients;                                                                                  % The estimate for the new data point
    end
elseif poly.Lim==1
%    yhat    = (poly.Coefficients'*[ones(1,size(X,2));X])';    
    yhat   = bsxfun(@plus,(poly.Coefficients(2:end,:)'*X),poly.Coefficients(1,:)')';
elseif poly.Lim==0
    yhat    = repmat(poly.Coefficients,size(X,2),1);
end

return