%% Demos to visualize uniform and adaptive partitions for 3D manifolds
%% First run ''RunExamples.m'', and then run ''Demo_PCMI.m'' for ExampleIdx==8,9,11,12,13,31

J        = max(gMRA.Scales); % maximal scale
LeafRow  = zeros(size(gMRA.cp)); for k=1:length(gMRA.LeafNodes); LeafRow(gMRA.LeafNodes(k)) = k; end
LenKappa = length(DataError_train.AdaptiveKappa);

%% display uniform approximations
figure; subplot(3,4,1); scatter3(X(1,:),X(2,:),X(3,:),'b'); axis tight; xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');  zLimits = get(gca,'ZLim');  %orginal data
if ExampleIdx>=11 && ExampleIdx<=13;  view(-25,-75); elseif ExampleIdx==31; view([0 90]); elseif ExampleIdx>=8 && ExampleIdx<=9; view([0 -1]); end
for j = 3:2:min(23,J)
    Projection  = GetApproximationOnPartition(gMRA,DataError_train.UniformPartition{j}.leafParentInPartition,XGWT);
    Label       = DataError_train.UniformPartition{j}.leafParentInPartition(LeafRow(gMRA.IniLabels));
    subplot(3,4,(j-1)/2+1); scatter3(Projection(1,:),Projection(2,:),Projection(3,:),[],Label); 
    axis([xLimits yLimits zLimits]); title(['j = ' num2str(j)]);
    if ExampleIdx>=11 && ExampleIdx<=13;  view(-25,-75); elseif ExampleIdx==31; view([0 90]); elseif ExampleIdx>=8 && ExampleIdx<=9; view([0 -1]); end
end
suptitle('GMRA')

%% display \|P_j x - P_{j+1} x\| at each point
[~, PointsOrder] = sort(TrainAndTestOpts.TrainIdxs,'ascend');
figure; imagesc(1:size(XGWT.DeltaPoint,1),1:length(DataError_train.UniformScales),log10((XGWT.DeltaPoint(PointsOrder,:)'))); colorbar
title('log_{10}||P_{j+1}x-P_j x||^2 '); xlabel('Point index');ylabel('Scale')

%% display adaptive approximation, color by partition
figure;         
subplot(3,4,1); scatter3(X(1,:),X(2,:),X(3,:),'b');
% xLimits = get(gca,'XLim'); yLimits = get(gca,'YLim');  zLimits = get(gca,'ZLim');  %orginal data
if ExampleIdx>=11 && ExampleIdx<=13;  view(-25,-75); elseif ExampleIdx==31; view([0 90]); end
for ikappa=1:3:min(36,LenKappa)
    Projection = GetApproximationOnPartition(gMRA,DataError_train.AdaptivePartition{ikappa}.leafParentInPartition,XGWT);
    Label       = DataError_train.AdaptivePartition{ikappa}.leafParentInPartition(LeafRow(gMRA.IniLabels));
    Scale       = DataError_train.AdaptivePartition{ikappa}.leafScale(LeafRow(gMRA.IniLabels));
    subplot(3,4,(ikappa-1)/3+1); scatter3(Projection(1,:),Projection(2,:),Projection(3,:),[],Label);
    axis([xLimits yLimits zLimits]); title(['\kappa = ' num2str(DataError_train.AdaptiveKappa(ikappa))]);
    if ExampleIdx>=11 && ExampleIdx<=13;  view(-25,-75); elseif ExampleIdx==31; view([0 90]); elseif ExampleIdx>=8 && ExampleIdx<=9; view([0 -1]); end
end
suptitle('Adaptive GMRA (color by partition)');



%% display adaptive approximation, color by scale
figure;cursplotidx=1;         
subplot(3,3,1); scatter3(X(1,:),X(2,:),X(3,:),'b'); if ExampleIdx>=11 && ExampleIdx<=13;  view(-25,-75); elseif ExampleIdx==31; view([0 90]); end
nlogradii     = round((-log10(DataError_train.CenterRadii))*100)/100;        ColorIndex    = jet(J);  PointsIndex   = cell(J,1);
for ikappa = 4:2:min(26,LenKappa)
    Projection = GetApproximationOnPartition(gMRA,DataError_train.AdaptivePartition{ikappa}.leafParentInPartition,XGWT);
    LenNodes    = length(DataError_train.AdaptivePartition{ikappa}.nodes);
    for k = 1 : LenNodes
        kscale                  = gMRA.Scales(DataError_train.AdaptivePartition{ikappa}.nodes(k));
        PointsIndex{kscale}     = [PointsIndex{kscale} XGWT.PointsInNet{DataError_train.AdaptivePartition{ikappa}.nodes(k)}];
    end
    for j = 1 : J
        if ~isempty(PointsIndex{j})
            colorindex = floor(J*(nlogradii(j)-nlogradii(1))/(nlogradii(end)-nlogradii(1)));
            if colorindex <= 0; colorindex = 1;  elseif colorindex > J; colorindex = J;  end
            subplot(3,3,cursplotidx); plot3(Projection(1,PointsIndex{j}),Projection(2,PointsIndex{j}),Projection(3,PointsIndex{j}),'o','MarkerSize',2,'color',ColorIndex(colorindex,:)); 
        end
        if j<J;  hold on;  end
    end
    if ExampleIdx>=11 && ExampleIdx<=13;  view(-25,-75); elseif ExampleIdx==31; view([0 90]); elseif ExampleIdx>=8 && ExampleIdx<=9; view([0 -1]); end
    axis([xLimits yLimits zLimits]);   title([' partition size = ' num2str(DataError_train.AdaptiveNCell(ikappa)) ]);
    caxis([nlogradii(1) nlogradii(end)]);  colormap(jet);
    cursplotidx = cursplotidx+1;
end
hp4 = get(subplot(3,3,3),'Position');
ct = colorbar('Ticks',linspace(nlogradii(1),nlogradii(end),J),'TickLabels',fliplr(nlogradii),...
            'Position', [hp4(1)+hp4(3)+0.02  0.1*hp4(2)  0.03  hp4(2)+hp4(3)*1]);
        set(get(ct,'title'),'string','-log_{10}radius');         
 suptitle('Adaptive GMRA (color by scale)');


% clearvars J Projection Label Scale LeafRow j ikappa LenKappa ct hp4 nlogradii colorindex ColorIndex PointsOrder
% clearvars k kscale PointsIndex
