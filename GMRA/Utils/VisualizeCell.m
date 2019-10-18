function VisualizeCell(gMRA,X,cellindex)
%% Visualize a cell
%% red: points in a cell; 
%% green: points not in the cell; 
%% black: center of the cell
%% Input: gMRA
%%        cellindex

J  = max(gMRA.Scales);       % maximal scale
j  = gMRA.Scales(cellindex); % scale of the cell

%%% center
center  = gMRA.Centers{cellindex};           % center of the cell

%% get points
%%% indices of the points in the cell
leafindex  = gMRA.CellIndex(gMRA.CellIndex(:,j)==cellindex,J); % leaf nodes of this cell
netPts     = find(ismember(gMRA.IniLabels,leafindex));     
%%% indices of the points outside the cell
netPtsRest = setdiff((1:size(X,2)),netPts); 

figure
scatter3(X(1,netPts),X(2,netPts),X(3,netPts),'ro')
hold on
scatter3(X(1,netPtsRest),X(2,netPtsRest),X(3,netPtsRest),'g*')
hold on
scatter3(center(1),center(2),center(3),'ko','filled')



