function DisplayRegressionErrorVsPartitionSize( X, gMRAReg, Y, Yhat, Yhat_adapt, unifError, adaptError, sigmaY, competitors )

Nfcns       = size(Y,2);
Ndegrees    = length(gMRAReg.RegressionOpts.degree);

%% Loop through the functions
% Create the figures
colors                  = {'r','b','k','m','c','g','r','b','k','m','c','g','r','b','k','m','c','g','r','b','k','m','c','g','r','b','k','m','c','g'};
degreecolors            = {'r','b','k','m','c','g','r','b','k','m','c','g','r','b','k','m','c','g','r','b','k','m','c','g','r','b','k','m','c','g'};
unifErrorLineStyle      = '-';
adaptErrorLineStyle     = '-';
unifErrorMarkerStyle    = 'o';
adaptErrorMarkerStyle   = 'd';
LineWidth               = 1;
MarkerSize              = 7;
ScatterMarkerSize       = 10;
ScatterColorMap         = 'hot';

lNsubplotsPerRow        = max([3,Ndegrees]);

for r = 1:Nfcns    
    figure
    MakeFigLargeAndNice( gcf );
    if isfield(gMRAReg.RegressionOpts,'FcnNames')
        try
            h=set(gcf,'Name',sprintf('Function %s',gMRAReg.RegressionOpts.FcnNames{r}));
            set(h,'interpreter','none');
        catch
            set(gcf,'Name',sprintf('Function %d',r));
        end
    else
        set(gcf,'Name',sprintf('Function %d',r));
    end
    subplot(3,lNsubplotsPerRow,1);
    colormap(ScatterColorMap);
    scatter3(X(1,:),X(2,:),X(3,:),ScatterMarkerSize,Y(:,r),'filled'); colorbar
    if isfield(gMRAReg.RegressionOpts,'FcnNames')
        try
            h=title(sprintf('Function %s',gMRAReg.RegressionOpts.FcnNames{r}));
            set(h,'interpreter','none');
        catch
            title(sprintf('Function number %d',r));
        end
    else
        title(sprintf('Function number %d',r));
    end
%    [minerror,minerroridx] = min(squeeze(unifError.rel(end,:,r)));
     error = unifError.rel(:,:,r);
     [minerror,minerroridx] = min(error(:));
     [minerror1,minerror2] = ind2sub(size(error),minerroridx);
    subplot(3,lNsubplotsPerRow,2);
    colormap(ScatterColorMap);
%    scatter3(X(1,:),X(2,:),X(3,:),ScatterMarkerSize,log10(abs(Yhat.Yhat{end}{minerroridx}(:,r)-Y(:,r))),'filled'); colorbar
%    title(sprintf('Error of best uniform approximation:\n degree %d, error=%f',gMRAReg.RegressionOpts.degree(end),minerror));
    scatter3(X(1,:),X(2,:),X(3,:),ScatterMarkerSize,log10(abs(Yhat.Yhat{minerror1}{minerror2}(:,r)-Y(:,r))),'filled'); colorbar
    title(sprintf('Error of best uniform approximation:\n degree %d, error=%f',gMRAReg.RegressionOpts.degree(minerror1),minerror));
%    [minerror,minerroridx] = min(squeeze(adaptError.rel(end,:,r)));
    error = adaptError.rel(:,:,r);
    [minerror,minerroridx] = min(error(:));
    [minerror1,minerror2] = ind2sub(size(error),minerroridx);
    subplot(3,lNsubplotsPerRow,3);
    colormap(ScatterColorMap);
%    scatter3(X(1,:),X(2,:),X(3,:),ScatterMarkerSize,log10(abs(Yhat_adapt.Yhat{end}{minerroridx}(:,r)-Y(:,r))),'filled');colorbar
%    title(sprintf('Error of best adaptive approximation:\n degree %d, error=%f',gMRAReg.RegressionOpts.degree(end),minerror));
    scatter3(X(1,:),X(2,:),X(3,:),ScatterMarkerSize,log10(abs(Yhat_adapt.Yhat{minerror1}{minerror2}(:,r)-Y(:,r))),'filled');colorbar
    title(sprintf('Error of best adaptive approximation:\n degree %d, error=%f',gMRAReg.RegressionOpts.degree(minerror1),minerror));
    % Display histograms for the scales of the adapative partitions
    subplot(3,lNsubplotsPerRow,3*lNsubplotsPerRow);
%    hist(gMRAReg.gMRA.Scales(Yhat_adapt.partition{end}{minerroridx,r}),1:max(gMRAReg.gMRA.Scales));
%    title(sprintf('Distribution of scales in adaptive partitions for decreasing kappa; degree %d',gMRAReg.RegressionOpts.degree(end)));
    hist(gMRAReg.gMRA.Scales(Yhat_adapt.partition{minerror1}{minerror2,r}),1:max(gMRAReg.gMRA.Scales));
    title(sprintf('Distribution of scales in adaptive partitions for decreasing kappa; degree %d',gMRAReg.RegressionOpts.degree(minerror1)));
    
    minPartitionSz = Inf; 
    maxPartitionSz = -Inf;    
    clear curaxes
    for i = 1:Ndegrees                                                                                                          % Display approximation error, on figure for every degree of local polynomials
        subplot(3,lNsubplotsPerRow,lNsubplotsPerRow+i);        
        if ~gMRAReg.RegressionOpts.isClassificationProblem
            yhatpartitionsz = cellfun('length',Yhat.partition{i});
            p=plot(yhatpartitionsz,log10(squeeze(unifError.rel(i,:,r))),colors{i});hold on;                                     % uniform approximation
            curaxes(i) = gca;
            set(p,'LineStyle',unifErrorLineStyle);
            set(p,'Marker',unifErrorMarkerStyle);
            set(p,'MarkerSize',MarkerSize);
            set(p,'LineWidth',LineWidth);
            set(p,'Color',degreecolors{i});
            yhatadaptpartitionsz = cellfun('length',Yhat_adapt.partition{i}(:,r));
            p=plot(yhatadaptpartitionsz,log10(squeeze(adaptError.rel(i,:,r))),colors{i});                                       % adaptive approximation
            set(p,'LineStyle',adaptErrorLineStyle);
            set(p,'Marker',adaptErrorMarkerStyle);
            set(p,'MarkerSize',MarkerSize);
            set(p,'LineWidth',LineWidth);
            set(p,'Color',degreecolors{i});
            legend({'Uniform','Adaptive'});axis tight;            
            title(sprintf('degree %d',gMRAReg.RegressionOpts.degree(i)));
            minPartitionSz = min([minPartitionSz,yhatpartitionsz]);
            maxPartitionSz = max([maxPartitionSz,yhatpartitionsz]);
            minPartitionSz = min([minPartitionSz;yhatadaptpartitionsz]);
            maxPartitionSz = max([maxPartitionSz;yhatadaptpartitionsz]);
            if sigmaY>0
                p=plot(minPartitionSz:maxPartitionSz,ones(maxPartitionSz-minPartitionSz+1,1)*log10(sigmaY));
                set(p,'LineStyle','--','LineWidth',2,'Color','k');
            end
        else
            %             figure
            %             plot(cellfun('length',Yhat.partition{r}),(squeeze(unifError.rel(i,:,r))),'b:o','LineWidth',2);hold on;
            %             plot(cellfun('length',Yhat_adapt.partition{r}),(squeeze(adaptError.rel(i,:,r))),'r-o','LineWidth',2);legend({'Uniform','Adaptive'});axis tight;
            %             title(sprintf('Classification error as a function of #cells for regression; degree %d, function %d',gMRAReg.RegressionOpts.degree(i),r));
        end
    end
    SyncXAxisOfSubplotsInFig( gcf, curaxes );
    SyncYAxisOfSubplotsInFig( gcf, curaxes );
        
    strings = {};
    clear curaxes
    subplot(3,lNsubplotsPerRow,2*lNsubplotsPerRow+1);    
    for i = 1:Ndegrees
%        if length(unique(Y))>length(Y)/10                                                                                      % MM: to fix in case of classification instead of regression...?
            p=plot(cellfun('length',Yhat.partition{i}),log10(squeeze(unifError.rel(i,:,r))),colors{i});hold on;                 % uniform approximation
            set(p,'LineStyle',unifErrorLineStyle);
            set(p,'Marker',unifErrorMarkerStyle);
            set(p,'MarkerSize',MarkerSize);
            set(p,'LineWidth',LineWidth);
            set(p,'Color',degreecolors{i});
            strings{i} = sprintf('degree %d',gMRAReg.RegressionOpts.degree(i));
%        end
    end
    
    for l = 1 : length(competitors.alg)
        p=plot(get(gca,'XTick'),ones(1,length(get(gca,'XTick')))*log10(competitors.err{l}(r)));
        set(p,'LineStyle','-','LineWidth',LineWidth,'Color',competitors.col{l});
        strings{l+Ndegrees} = sprintf(competitors.alg{l},competitors.err{l});
        text()
    end
    
    if sigmaY>0
        p=plot(get(gca,'XTick'),ones(1,length(get(gca,'XTick')))*log10(sigmaY));
        set(p,'LineStyle','--','LineWidth',LineWidth,'Color','k');
    end
    title('Uniform partitions');
    legend(strings);
    axis tight;    
    curaxes(1) = gca;
   
    clear strings
    subplot(3,lNsubplotsPerRow,2*lNsubplotsPerRow+2);
    for i = 1:Ndegrees
        if ~gMRAReg.RegressionOpts.isClassificationProblem
            p=plot(cellfun('length',Yhat_adapt.partition{i}(:,r)),log10(adaptError.rel(i,:,r)));hold on   
            set(p,'LineStyle',adaptErrorLineStyle);
            set(p,'Marker',adaptErrorMarkerStyle);
            set(p,'MarkerSize',MarkerSize);
            set(p,'LineWidth',LineWidth);
            set(p,'Color',degreecolors{i});
            strings{i} = sprintf('degree %d',gMRAReg.RegressionOpts.degree(i));           
        end
    end
%     for l = 1 : length(competitors.err)
%         p=plot(get(gca,'XTick'),ones(1,length(get(gca,'XTick')))*log10(competitors.err{l}));
%         set(p,'LineStyle','-','LineWidth',LineWidth,'Color',competitors.col{l});
%     end
    
    for l = 1 : length(competitors.alg)
        p=plot(get(gca,'XTick'),ones(1,length(get(gca,'XTick')))*log10(competitors.err{l}(r)));
        set(p,'LineStyle','-','LineWidth',LineWidth,'Color',competitors.col{l});
        strings{l+Ndegrees} = sprintf(competitors.alg{l},competitors.err{l});
        text()
    end
    title('Adaptive partition');
    if sigmaY>0
        p=plot(get(gca,'XTick'),ones(1,length(get(gca,'XTick')))*log10(sigmaY));
        set(p,'LineStyle','--','LineWidth',LineWidth,'Color','k');
    end    
    legend(strings);
    axis tight;
    curaxes(2) = gca;    
    
    SyncXAxisOfSubplotsInFig( gcf, curaxes );
    SyncYAxisOfSubplotsInFig( gcf, curaxes );
    
    
%    clear strings            
%     figure;
%     for i = 1:length(gMRAReg.RegressionOpts.degree)
%         if length(unique(Y))>length(Y)/10
%             plot(cellfun('length',Yhat_adapt.partition{r}),log10(squeeze(adaptError.rel(i,:,r))),'r-o','LineWidth',LineWidth);axis tight;hold on;
%             title(sprintf('Relative error as a function of #cells for regression, adaptive, function %d',r));
%         else
%             plot(cellfun('length',Yhat_adapt.partition{i}(:,r)),adaptError.rel(i,:,r),'r-o','LineWidth',2);legend({'Uniform','Adaptive'});axis tight;title('Classification error as a function of #cells for regression');
%         end
%     end
end


return
