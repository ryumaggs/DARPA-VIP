function AGWT_DisplayApproxErr(gMRA,DataError,displayopts)
%
% Function to display approximation error versus scale/partition size 
%
%
%
% Input: gMRA:         GMRA data structure
%        DataError:    structure containing errors of uniform and adaptive GMRA
%        displayopts:  data structure of options:
%                   [errortype]: 0 absolute error  1 relative error
%                   [errornorm]: '2' L^2 error 'inf' L^\infty error
%                   [noisesize]: plot an horizontal line at this level if provided and >0
%
% (c) Wenjing Liao, Mauro Maggioni

if nargin < 3
    displayopts = struct('errortype',0,'errornorm',2);                     % default: absolute error; L2 norm
end
if ~isfield(displayopts,'noisesize'),   displayopts.noisesize = 0;  end;

ErrOpts = displayopts.errortype;
switch displayopts.errornorm
    case '2';        ErrCol = 1;
    case 'inf';      ErrCol = 2;
end


J  = max(gMRA.Scales);  % maximum scale
d  = gMRA.opts.ManifoldDimension; % intrinsic dimension

%% Determine the scales to fit in GMRA
if ErrOpts == 0
    goodidxs                    = intersect(find(isfinite(DataError.UniformScales)),find(isfinite(DataError.UniformAbsolute(:,ErrCol))));
    goodidxs([1,2,end-1,end])   = [];
    goodidxs                    = intersect(goodidxs,find(diff(DataError.UniformAbsolute(:,ErrCol))<median(diff(DataError.UniformAbsolute(:,ErrCol)))/2));
elseif ErrOpts == 1
    goodidxs                    = intersect(find(isfinite(DataError.UniformScales)),find(isfinite(DataError.UniformRelative(:,ErrCol))));
    goodidxs([1,2,end-1,end])   = [];
    goodidxs                    = intersect(goodidxs,find(diff(DataError.UniformRelative(:,ErrCol))<median(diff(DataError.UniformRelative(:,ErrCol)))/2));
    %goodidxs = find(DataError.UniformAbsolute(:,ErrCol)>2*eps*max(DataError.UniformAbsolute(:,ErrCol)));
end
UniformFitRange = goodidxs(1):goodidxs(end);
% fprintf('Uniform  ls fit range: %6.0f to %6.0f out of 1 to %6.0f\n',goodidxs(1),goodidxs(end),J)


%% Error versus scale - uniform GMRA
if ~isfield(gMRA,'CoverTree') || ~isfield(gMRA.CoverTree,'theta'), gMRA.CoverTree.theta = 2;    end
if ErrOpts == 0
    % Absolute Error: Least square fitting    
    GMRAUniformFit    = polyfit(-log10(DataError.UniformRadii(UniformFitRange)),log10(DataError.UniformAbsolute(UniformFitRange,ErrCol)),1);
    CenterUniformFit  = polyfit(-log10(DataError.CenterRadii(UniformFitRange)),log10(DataError.CenterAbsolute(UniformFitRange,ErrCol)),1);
    figure
    set(gcf,'Name','Absolute error');
    plot(-log10(DataError.CenterRadii),log10(DataError.CenterAbsolute(:,ErrCol)),'g*-.');    hold on;
    pc = plot(-log10(DataError.CenterRadii(UniformFitRange)),log10(DataError.CenterAbsolute(UniformFitRange,ErrCol)),'g*-');    hold on;
    plot(-log10(DataError.UniformRadii),log10(DataError.UniformAbsolute(:,ErrCol)),'bo-.'); hold on;
    pg = plot(-log10(DataError.UniformRadii(UniformFitRange)),log10(DataError.UniformAbsolute(UniformFitRange,ErrCol)),'bo-');
    % mark percentage of cells whose size >= d
    %for j = 1 : J
    %    if gMRA.opts.ManifoldDimension > 0
    %        text(double(DataError.UniformScales(j)),double(log10(DataError.UniformAbsolute(j,1))/log10(1/gMRA.CoverTree.theta)),...
    %         num2str(100*DataError.PercentBd(j),'%6.0f'))
    %    end
    %    text(double(DataError.CenterScales(j)),double(log10(DataError.CenterAbsolute(j,1))/log10(1/gMRA.CoverTree.theta)),...
    %         num2str(j,'%6.0f'))
    %end
    legend([pc pg],['Center: slope= ',num2str(CenterUniformFit(1))],['GMRA:  slope= ',num2str(GMRAUniformFit(1))])
    ylabel('log_{10}(error)')
elseif ErrOpts == 1
    % Relative Error: Least square fitting    
    GMRAUniformFit    = polyfit(-log10(DataError.UniformRadii(UniformFitRange)),log10(DataError.UniformRelative(UniformFitRange,ErrCol)),1);
    CenterUniformFit  = polyfit(-log10(DataError.CenterRadii(UniformFitRange)),log10(DataError.CenterRelative(UniformFitRange,ErrCol)),1);
    figure
    set(gcf,'Name','Relative error');
    plot(-log10(DataError.CenterRadii),log10(DataError.CenterRelative(:,ErrCol)),'g*-.');    hold on;
    pc = plot(-log10(DataError.CenterRadii(UniformFitRange)),log10(DataError.CenterRelative(UniformFitRange,ErrCol)),'g*-');    hold on;
    plot(-log10(DataError.UniformRadii),log10(DataError.UniformRelative(:,ErrCol)),'bo-.'); hold on;
    pg = plot(-log10(DataError.UniformRadii(UniformFitRange)),log10(DataError.UniformRelative(UniformFitRange,ErrCol)),'bo-');
    legend([pc pg],['Center: slope= ',num2str(CenterUniformFit(1))],['GMRA:  slope= ',num2str(GMRAUniformFit(1))])
    ylabel('log_{10}(relative error)')    
end
if displayopts.noisesize>0
    plot(-log10(DataError.UniformRadii),log10(displayopts.noisesize*ones(size(DataError.UniformRadii))),'k--'); hold on;
end
title('log_{10}(error) versus -log_{10}(radii)')
xlabel('-log_{10}(radii)')
grid on
% set(gca,'FontSize',14,'FontWeight','bold')


%% Determine the range to fit in adaptive GMRA
AdaptiveFitRange                    = find((DataError.AdaptiveNCell>=DataError.UniformNCell(UniformFitRange(1))).*(DataError.AdaptiveNCell<=DataError.UniformNCell(UniformFitRange(end))));
AdaptiveFitRange([1,2,end-1,end])   = [];
% fprintf('Adaptive ls fit range: %6.0f to %6.0f out of 1 to %6.0f\n',AdaptiveFitRange(1),AdaptiveFitRange(end),length(DataError.AdaptiveKappa))


%% Error versus partition size - uniform + adaptive
figure
if ErrOpts == 0
    UniformFit    = polyfit(log10(DataError.UniformNCell(UniformFitRange))',log10(DataError.UniformAbsolute(UniformFitRange,ErrCol)),1);
    AdaptiveFit   = polyfit(log10(DataError.AdaptiveNCell(AdaptiveFitRange)),log10(DataError.AdaptiveAbsolute(AdaptiveFitRange,ErrCol)),1);
    plot(log10(DataError.UniformNCell)',log10(DataError.UniformAbsolute(:,ErrCol)),'bo-.'); hold on;
    pu = plot(log10(DataError.UniformNCell(UniformFitRange))',log10(DataError.UniformAbsolute(UniformFitRange,ErrCol)),'bo-'); hold on;
    plot(log10(DataError.AdaptiveNCell),log10(DataError.AdaptiveAbsolute(:,ErrCol)),'ro-.'); hold on;
    pa = plot(log10(DataError.AdaptiveNCell(UniformFitRange)),log10(DataError.AdaptiveAbsolute(UniformFitRange,ErrCol)),'ro-');
    legend([pu pa],['Uniform:  slope= ',num2str(UniformFit(1))],...
           ['Adaptive: slope= ',num2str(AdaptiveFit(1))])
    ylabel('log_{10}(error)')
elseif ErrOpts == 1   
  UniformFit    = polyfit(log10(DataError.UniformNCell(UniformFitRange))',log10(DataError.UniformRelative(UniformFitRange,ErrCol)),1);
    AdaptiveFit   = polyfit(log10(DataError.AdaptiveNCell(AdaptiveFitRange)),log10(DataError.AdaptiveRelative(AdaptiveFitRange,ErrCol)),1);
    plot(log10(DataError.UniformNCell)',log10(DataError.UniformRelative(:,ErrCol)),'bo-.'); hold on;
    pu = plot(log10(DataError.UniformNCell(UniformFitRange))',log10(DataError.UniformRelative(UniformFitRange,ErrCol)),'bo-'); hold on;
    plot(log10(DataError.AdaptiveNCell),log10(DataError.AdaptiveRelative(:,ErrCol)),'ro-.'); hold on;
    pa = plot(log10(DataError.AdaptiveNCell(AdaptiveFitRange)),log10(DataError.AdaptiveRelative(AdaptiveFitRange,ErrCol)),'ro-');
    legend([pu pa],['Uniform:  slope= ',num2str(UniformFit(1))],...
           ['Adaptive: slope= ',num2str(AdaptiveFit(1))])
    ylabel('log_{10}(relative error)')
end
if displayopts.noisesize>0
    plot(log10(DataError.UniformNCell),log10(displayopts.noisesize*ones(size(DataError.UniformNCell))),'k--'); hold on;
end
title('log_{10}(error) versus log_{10}(partition size)')
xlabel('log_{10}(partition size)')
%set(gca,'slope',14,'FontWeight','bold')
grid on 


% %% Error versus entropy - uniform + adaptive
% figure
% if ErrOpts == 0
%     UniformFit    = polyfit(log10(DataError.UniformPartitionEntropy(UniformFitRange)),log10(DataError.UniformAbsolute(UniformFitRange,ErrCol)),1);
%     AdaptiveFit   = polyfit(log10(DataError.AdaptivePartitionEntropy(AdaptiveFitRange)),log10(DataError.AdaptiveAbsolute(AdaptiveFitRange,1)),1);
%     plot(log10(DataError.UniformPartitionEntropy)',log10(DataError.UniformAbsolute(:,ErrCol)),'bo--','MarkerSize',4,'MarkerFaceColor','b'); hold on;
%     plot(log10(DataError.AdaptivePartitionEntropy),log10(DataError.AdaptiveAbsolute(:,ErrCol)),'ro-','MarkerSize',4,'MarkerFaceColor','r')
%     legend(['Uniform:  slope= ',num2str(UniformFit(1))],...
%            ['Adaptive: slope= ',num2str(AdaptiveFit(1))])
%     ylabel('log_{2}(error)')
% elseif ErrOpts == 1   
%     UniformFit    = polyfit(log10(DataError.UniformPartitionEntropy(UniformFitRange)),log10(DataError.UniformRelative(UniformFitRange,ErrCol)),1);
%     AdaptiveFit   = polyfit(log10(DataError.AdaptivePartitionEntropy(AdaptiveFitRange)),log10(DataError.AdaptiveRelative(AdaptiveFitRange,ErrCol)),1);
%     plot(log10(DataError.UniformPartitionEntropy)',log10(DataError.UniformRelative(:,ErrCol)),'bo--','MarkerSize',4,'MarkerFaceColor','b'); hold on;
%     plot(log10(DataError.AdaptivePartitionEntropy),log10(DataError.AdaptiveRelative(:,ErrCol)),'ro-','MarkerSize',4,'MarkerFaceColor','r')
%     legend(['Uniform:  slope= ',num2str(UniformFit(1))],...
%            ['Adaptive: slope= ',num2str(AdaptiveFit(1))])
%     ylabel('log_{2}(relative error)')
% end
% title('log_{2}(error) versus log_{2}(entropy)')
% xlabel('log_{2}(entropy)')
% %set(gca,'slope',14,'FontWeight','bold')
% grid on 
% 
% 
% 
% %% Error versus entropy - adaptive
% figure
% if ErrOpts == 0
%     if d >= 3
%         LogEntropy              = log10(DataError.AdaptivePartitionEntropy);
%         LogErr                  = (d-2)*log10(DataError.AdaptiveAbsolute(:,ErrCol));
%         AdaptiveFitErrEntropy   = polyfit(LogEntropy,LogErr,1);
%         LogEta                  = (d-2)*log10(DataError.AdaptiveKappa)+(d-2)/2*log10(DataError.AdaptivePartitionEntropy);
%         AdaptiveFitEtaEntropy   = polyfit(LogEntropy,LogEta,1);
%         plot(LogEntropy,LogEta,'bo-',LogEntropy,LogErr,'r*-','MarkerSize',4)
%         legend(['\eta \sim entropy:       slope= ',num2str(AdaptiveFitEtaEntropy(1))],...
%            ['error \sim entropy: slope= ',num2str(AdaptiveFitErrEntropy(1))])
%     else
%         LogEntropy              = log10(DataError.AdaptivePartitionEntropy);
%         LogErr                  = log10(DataError.AdaptiveAbsolute(:,ErrCol));
%         AdaptiveFitErrEntropy   = polyfit(LogEntropy,LogErr,1);
%         plot(LogEntropy,LogErr,'r*-','MarkerSize',4)
%         legend(['error \sim entropy: slope= ',num2str(AdaptiveFitErrEntropy(1))])
%     end
% elseif ErrOpts == 1
%     if d >= 3
%         LogEntropy              = log10(DataError.AdaptivePartitionEntropy);
%         LogErr                  = (d-2)*log10(DataError.AdaptiveRelative(:,ErrCol));
%         AdaptiveFitErrEntropy   = polyfit(LogEntropy,LogErr,1);
%         LogEta                  = (d-2)*log10(DataError.AdaptiveKappa)+(d-2)/2*log10(DataError.AdaptivePartitionEntropy);
%         AdaptiveFitEtaEntropy   = polyfit(LogEntropy,LogEta,1);
%         plot(LogEntropy,LogEta,'bo-',LogEntropy,LogErr,'r*-','MarkerSize',4)
%         legend(['\eta \sim entropy:       slope= ',num2str(AdaptiveFitEtaEntropy(1))],...
%            ['relative error \sim entropy: slope= ',num2str(AdaptiveFitErrEntropy(1))])
%     else
%         LogEntropy              = log10(DataError.AdaptivePartitionEntropy);
%         LogErr                  = log10(DataError.AdaptiveRelative(:,ErrCol));
%         AdaptiveFitErrEntropy   = polyfit(LogEntropy,LogErr,1);
%         plot(LogEntropy,LogErr,'r*-','MarkerSize',4)
%         legend(['relative error \sim entropy: slope= ',num2str(AdaptiveFitErrEntropy(1))])
%     end
% end
% if d >= 3
%     title('log_{10}(\eta^{d-2}entropy^{(d-2)/2}) and log_{10}(error^{d-2}) versus log_{10}(entropy)')
% else
%     title('log_{10}(error) versus log_{10}(entropy)')
% end
% xlabel('log_{10}(entropy)')
% % set(gca,'FontSize',14,'FontWeight','bold')
% grid on 
% 
% 

