function [handle,fit]=PlotPolyFit( X, Y, N, Sigma )

if nargin<3 || isempty(N),      N=1;           end
if nargin<4 || isempty(Sigma),  Sigma=0;       end

goodidxs            = intersect(find(isfinite(X)),find(isfinite(Y)));
goodidxs([1,end])   = [];

if Sigma~=0
    goodidxs = intersect(goodidxs,find(Y>Sigma));
end

new_goodidxs    = intersect(goodidxs,find(diff(Y)<median(diff(Y))/2));
if length(new_goodidxs)>2
    goodidxs = new_goodidxs;
end
badidxs     = setdiff(1:length(X),goodidxs);
fit         = polyfit(X(goodidxs),Y(goodidxs),N);

hold on;
handle      = plot(X,polyval(fit,X),'g--', 'LineWidth', 1);
plot(X(badidxs),Y(badidxs),'x');

legend(handle,['y=',poly2str(fit,'x')],'Location','Best');

return