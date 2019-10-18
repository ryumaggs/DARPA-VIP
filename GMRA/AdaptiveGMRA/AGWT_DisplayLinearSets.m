function AGWT_DisplayLinearSets(gMRA,X,XGWT,DataError)

J = max(gMRA.Scales);

%% Plot X
figure
scatter3(X(1,:),X(2,:),X(3,:),'b')
axis tight
xLimits = get(gca,'XLim');  %# Get the range of the x axis
yLimits = get(gca,'YLim');  %# Get the range of the y axis
zLimits = get(gca,'ZLim');  %# Get the range of the z axis
view(30,70);


% Plot linear approximations
figure
for j = 1:J
    Projectionj = GetApproximationOnPartition(gMRA,DataError.UniformPartition{j}.leafParentInPartition,XGWT);
    scatter3(Projectionj(1,:),Projectionj(2,:),Projectionj(3,:),'b')
    axis([xLimits yLimits zLimits])
    title(['j = ' num2str(j)])
    view(30,70);
    pause;
end



% j=4;
% Projectionj = IGWTScalej(gMRA,XGWT,j);
% subplot(1,3,2)
% scatter3(Projectionj(1,:),Projectionj(2,:),Projectionj(3,:),'b')
% axis([xLimits yLimits zLimits])
% title(['j = ' num2str(j)])
% view(16,74);
% 
% j=8;
% Projectionj = IGWTScalej(gMRA,XGWT,j);
% subplot(1,3,3)
% scatter3(Projectionj(1,:),Projectionj(2,:),Projectionj(3,:),'b')
% axis([xLimits yLimits zLimits])
% title(['j = ' num2str(j)])
% view(16,74);


% 
% %% Adaptive GMRA
% AdaptiveAbsError = sqrt(norm((X - DataAGWT.Projection).^2,'fro')/size(X,2));
% AdaptiveReError  = AdaptiveAbsError/(sqrt(norm(X.^2,'fro')/size(X,2)));
% AdaptiveNCell = length(unique(DataAGWT.PartitionIndex));
% fprintf('Adaptive: kappa=%6.2f Cell=%6.0f AbsError=%e ReError=%6.2f\n',DataAGWT.kappa,AdaptiveNCell,AdaptiveAbsError,AdaptiveReError)
% 
% 
% figure
% scatter3(DataAGWT.Projection(1,:),DataAGWT.Projection(2,:),DataAGWT.Projection(3,:))
% title(['kappa = ' num2str(DataAGWT.kappa) ' N = ' num2str(AdaptiveNCell) ' Error = ' num2str(AdaptiveAbsError) ])
% axis tight
% view([0 10]);
% 
% %% Plot Deltajk at each point
% %figure
% %imagesc(log2(DataAGWT.DeltaPoint.'))
% %colorbar
% %set(gca, 'CLim', [-5 1]);
% 


% %Plot points at each cell
% figure
% for i = 1: nLeafNodes
%   iFineNet = gMRA.LeafNodes(i);
%   j = gMRA.Scales(iFineNet);
%   netPtsidx = find(gMRA.IniLabels == iFineNet);
%   netPts    = X(:,netPtsidx);
%   clf
%   scatter3(X(1,:),X(2,:),X(3,:),'g')
%   hold on
%   scatter3(netPts(1,:),netPts(2,:),netPts(3,:),'k','filled','d')
%   axis([xLimits yLimits zLimits])
%   title(num2str(i))
%   view([0 10]);
%   CellM(i)=getframe(gcf);
%   %% pause
% end
% 
% movie2avi(CellM,'CellMovie.avi');



% % 
% %% Plot points at each cell
% figure
% for i = 1: 200
%     index = i*100:(i+1)*100;
%     scatter3(X(1,index),X(2,index),X(3,index))
%     axis([xLimits yLimits zLimits])
%     view([16 74])
%     pause
% end