function AGWT_DisplayBallRadii(gMRA,BallInCell,d)
%% Input: BallInCell -- The approximate Radii of the largest ball contained in every cell 
%% The statistics of BallInCell is studied.

J = length(BallInCell);
% Radii in cover tree
CoverTreeRadii = gMRA.CoverTree.radii;
Jc             = min(length(CoverTreeRadii),J);


% Mean
MeanRadii         = zeros(1,J);  % mean of in radii
MedianRadii       = zeros(1,J);  % median of in radii
UniformRadii      = zeros(J,1);  % mean of outer radii
StdRadii          = zeros(1,J);
QuantileRadii95   = zeros(1,J);
QuantileRadii05   = zeros(1,J);
RatioOuterInRadii = BallInCell;
MeanRatio         = zeros(1,J);
StdRatio          = zeros(1,J);
Scales            = -log2(CoverTreeRadii); %1:J;


for k = 1:J
%     if nargin == 3
%         BallInCell{k}(:,4) = BallInCell{k}(:,2).^d./BallInCell{k}(:,3);
%     end
    MeanRadii(k)        = mean(BallInCell{k}(:,2));
    MedianRadii(k)      = median(BallInCell{k}(:,2));
    StdRadii(k)         = std(BallInCell{k}(:,2));
    QuantileRadii95(k)  = quantile(BallInCell{k}(:,2),0.95);
    QuantileRadii05(k)  = quantile(BallInCell{k}(:,2),0.05);
    UniformRadii(k)     = mean(gMRA.Radii(gMRA.Scales==k));%mean(gMRA.Radii(uniquev));
    Nnode               = size(BallInCell{k},1);
    for l = 1:Nnode
        RatioOuterInRadii{k}(l,2) = gMRA.Radii(RatioOuterInRadii{k}(l,1))/BallInCell{k}(l,2);
    end
    MeanRatio(k)        =  sum(RatioOuterInRadii{k}(:,2).*RatioOuterInRadii{k}(:,3));
    StdRatio(k)         =  std(RatioOuterInRadii{k}(:,2));
end
FitRange       = 1:J;
OuterFit       = polyfit(Scales(FitRange),log2((UniformRadii(FitRange)))',1);
MeanFit        = polyfit(Scales(FitRange),log2(MeanRadii(FitRange)),1);
MedianFit      = polyfit(Scales(FitRange),log2(MedianRadii(FitRange)),1);
Quantile95Fit  = polyfit(Scales(FitRange),log2(QuantileRadii95(FitRange)),1);
Quantile05Fit  = polyfit(Scales(FitRange),log2(QuantileRadii05(FitRange)),1);

if isempty(FitRange)
    fprintf('Not enough scales to compute the in-radii\n')
    return;
end


% figure
% plot(Scales(FitRange),log2(UniformRadii(FitRange)),'k+-',Scales(FitRange),log2(QuantileRadii95(FitRange)),'gp-',Scales(FitRange),log2(MeanRadii(FitRange)),'r*-',Scales(FitRange),log2(MedianRadii(FitRange)),'bo-',Scales(FitRange),log2(QuantileRadii05(FitRange)),'md-')
% legend(['outer-radii:      slope = ',num2str(OuterFit(1))],['in-radii 0.95 quantile: slope = ',num2str(Quantile95Fit(1))],['in-radii mean:             slope = ',num2str(MeanFit(1))],['in-radii median:          slope = ',num2str(MeanFit(1))],['in-radii 0.05 quantile: slope = ',num2str(Quantile05Fit(1))])
% xlabel('scale')
% ylabel('log2(radius)')
% title('log2(radii) versus scale')

figure
plot(-log2(CoverTreeRadii(FitRange)),log2(UniformRadii(FitRange)),'k+-',-log2(CoverTreeRadii(FitRange)),log2(QuantileRadii95(FitRange)),'gp-',-log2(CoverTreeRadii(FitRange)),log2(MeanRadii(FitRange)),'r*-',-log2(CoverTreeRadii(FitRange)),log2(MedianRadii(FitRange)),'bo-',-log2(CoverTreeRadii(FitRange)),log2(QuantileRadii05(FitRange)),'md-')
legend(['outer-radii:      slope = ',num2str(OuterFit(1))],['in-radii 0.95 quantile: slope = ',num2str(Quantile95Fit(1))],['in-radii mean:             slope = ',num2str(MeanFit(1))],['in-radii median:          slope = ',num2str(MeanFit(1))],['in-radii 0.05 quantile: slope = ',num2str(Quantile05Fit(1))])
xlabel('scale')
ylabel('log2(radius)')
title('log2(radii) versus scale j')
ylim([-5 1])
xlim([-log2(CoverTreeRadii(FitRange(1))) -log2(CoverTreeRadii(FitRange(end)))])

figure;plot(1:Jc,log2(CoverTreeRadii(1:Jc)),'r*--',1:Jc,log2(UniformRadii(1:Jc)),'bo-',1:Jc,log2(MeanRadii(1:Jc)),'gp-.')
legend('Cover tree radii','Mean of outer-radii','Mean of in-radii')
title('log2(radii) versus scale')
xlabel('scale')
ylabel('log2(radii)')

figure
errorbar(-log2(CoverTreeRadii(FitRange)),MeanRatio(FitRange),StdRatio(FitRange),'bp')
title('Ratio between outer-radii and in-radii')
xlabel('scale')





% if nargin == 4
%     figure
%     for k = 1:J
%         plot(k,log10(BallInCell{k}(:,4)),'*','MarkerSize',8); hold on;
%     end
%     title('Ratio between r^d and empirical measure for all cells')
%     xlabel('scale')
%     ylabel('ratio')
%     set(gca,'fontweight','bold','fontsize',14)
% end
% 
% figure
% for k = 1:J
%     plot(k,log10(BallInCell{k}(:,4)/(pi^2*gMRA.TreeOnLeaves.radii(k)^gMRA.opts.ManifoldDimension)),'*','MarkerSize',8); hold on;
% end

% figure;plot(MeanRadii/(pi^2))
% figure;plot(MedianRadii/(pi^2))
% figure;plot(StdRadii/(pi^2))
% figure;plot(StdRadii./MeanRadii)
%figure;plot(QuantileRadii/(pi^2))
% 
% figure;loglog(gMRA.TreeOnLeaves.radii(1:J),MeanRadii);
% 
% figure;plot(MeanRadii./(gMRA.TreeOnLeaves.radii(1:J)));
% figure;plot(MedianRadii./(gMRA.TreeOnLeaves.radii(1:J)));


%figure;plot(QuantileRadii./MeanRadii);
% for k = 1:length(BallInCell);meanRadii(k)=mean(BallInCell{k}(:,2).^2./BallInCell{k}(:,3));medianRadii(k)=median(BallInCell{k}(:,2).^2./BallInCell{k}(:,3));stdRadii(k)=std(BallInCell{k}(:,2).^2./BallInCell{k}(:,3));quantileRadii(k)=quantile(BallInCell{k}(:,2).^2./BallInCell{k}(:,3),0.05);end
% figure;plot(meanRadii)
% figure;plot(meanRadii/(pi^2))
% figure;plot(medianRadii/(pi^2))
% figure;plot(stdRadii/(pi^2))
% figure;plot(stdRadii./meanRadii)
% figure;plot(quantileRadii/(pi^2))
% help quantile
% help prctile
% prctile(rand(100,1),0.9)
% prctile(rand(100,1),90)
% help quantile
% for k = 1:length(BallInCell);meanRadii(k)=mean(BallInCell{k}(:,2).^2./BallInCell{k}(:,3));medianRadii(k)=median(BallInCell{k}(:,2).^2./BallInCell{k}(:,3));stdRadii(k)=std(BallInCell{k}(:,2).^2./BallInCell{k}(:,3));quantileRadii(k)=quantile(BallInCell{k}(:,2).^2./BallInCell{k}(:,3),0.9);end
% figure;plot(quantileRadii/(pi^2))
% for k = 1:length(BallInCell);meanRadii(k)=mean(BallInCell{k}(:,2).^2./BallInCell{k}(:,3));medianRadii(k)=median(BallInCell{k}(:,2).^2./BallInCell{k}(:,3));stdRadii(k)=std(BallInCell{k}(:,2).^2./BallInCell{k}(:,3));quantileRadii(k)=quantile(BallInCell{k}(:,2).^2./BallInCell{k}(:,3),0.99);end
% figure;plot(quantileRadii/(pi^2))