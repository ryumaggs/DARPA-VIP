function [figHandle,plotHandles,YLim,fit] = DisplayApproximationErrors( gMRA, ApproxError, Quantiles, Sigma, figHandle, Color, LineStyle, theta )

if nargin<3 Quantiles = [];                             end
if nargin<4 || isempty(Sigma),      Sigma = 0.0;        end
if nargin<5 || isempty(figHandle),  figHandle=figure;   end
if nargin<6 || isempty(Color),      Color = [0,0,0];    end
if nargin<7 || isempty(LineStyle),  LineStyle = '--';    end
if nargin<8 || isempty(theta),      theta = gMRA.TreeOnLeaves.theta;  end

radii = gMRA.Radii_AvePerScale(1:size(ApproxError,1));

figure(figHandle);
[plotHandles(1),fit] = PlotPolyFit ( log10(1./radii')/log10(1/theta), log10(ApproxError(:,1))/log10(1/theta),[],...             % Plot fit to approximation error
                                     log10(2*Sigma)/log10(1/theta));
set(plotHandles(1),'Color',Color*0.9,'LineWidth',0.5,'LineStyle',LineStyle);

if ~isempty(Quantiles)                                                                                                          % Plot quantiles if provided
    hold on;
    colors = linspace(0.75,0.25,size(Quantiles,3));
    for l = 1:size(Quantiles,3)
        plot(log10(1./radii)/log10(1/theta),log10(Quantiles(:,1,l))/log10(1/theta),LineStyle,'Color',(1-Color)*colors(l));
    end
end

plotHandles(2) = plot( log10(1./radii)/log10(1/theta), log10(ApproxError(:,1))/log10(1/theta),'o-','LineWidth',1,'Color',Color); % Plot approximation error

xlabel('log_{1/\theta}(1/r)');ylabel('log_{1/\theta}(Relative L^2 error)');title('Relative approximation error in L^2');axis tight;grid on;

curYLim = get(gca,'YLim');
if ~isempty(Quantiles),
    set(gca,'YLim',[min([min(min((log10(Quantiles(:,1,:))/log10(1/theta)))),min(log10(ApproxError(:,1))/log10(1/theta))]),curYLim(2)]);
else
    set(gca,'YLim',[min(log10(ApproxError(:,1))/log10(1/theta)),max(log10(ApproxError(:,1))/log10(1/theta))]);
end
YLim = get(gca,'YLim');

if Sigma>0                                                                                                                      % Display noise level
    plot(log10(1./radii)/log10(1/theta),ones(1,length(radii))*log10(Sigma)/log10(1/theta),'Color',[0.75,0.1,0.1],'LineStyle',':');
    yticks=[get(gca,'YTick'),round(100*log10(Sigma)/log10(1/theta))/100];
    set(gca,'YTick',unique(yticks));
end


return
