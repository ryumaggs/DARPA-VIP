function BallInCell = ComputeRadiiBallInCell(gMRA,X,tol,Times)
%% Approximately compute the radius of the largst ball that every cell contains
%%% Input: gMRA: GMRA structure
%%%        X:    data used to construct GMRA
%%%        tol:  percentage of points required in a cell
%%%        Times:number of centers to search the ball
%%% For example: BallInCell = ComputeRadiiBallInCell(gMRA,X,0.99,5)
%%% computes the largest radius of the ball in every cell such that 99% 
%%% points in the ball is from the cell, and we search the ball for 5
%%% centers.
%%%
%%% Output: BallInCell is a structure of length J(the maximal scale)
%%%         BallInCell{j} has three columns.
%%%                       First column: Index of the cell
%%%                       Second column: In-radius of the ball
%%%                       Third column: Empirical measure of the cell

if nargin ==3 
    Times = 1;
elseif nargin == 2
    Times = 1;
    tol   = 0.95;
end

J          = max(gMRA.Scales);
nLeafNodes = numel(gMRA.LeafNodes); % number of leaf nodes
nAllPts    = size(X,2);

CellIndex   = zeros(nLeafNodes,J);  %cell index      % cell index
%% Label cell index to easily search the corresponding leaf nodes of any cell
for i = 1: nLeafNodes
    iFineNet = gMRA.LeafNodes(i);    
    j = gMRA.Scales(iFineNet);
    while j>=1
        CellIndex(i,j) = iFineNet;  
        % pass to the next scale
        j = j-1;
        iFineNet   = gMRA.cp(iFineNet);
    end
end

J          = find(min(CellIndex,[],1),1,'last');
BallInCell = cell(1,J);



for j = 1:J
    CellIndexj = unique(CellIndex(:,j));
    NCellj     = length(CellIndexj);
    BallInCell{j} = zeros(NCellj,3);           % First column: cell index; 
                                               % Second: radii; 
                                               % Third: empirical measure of the cell
    fprintf('Checking cells at scale j=%6.0f \n',j)
    for iter = 1 : NCellj
        cellindex    = CellIndexj(iter);                     % cell index
        %% get points
        %%% indices of the points in the cell
        leafindex    = CellIndex(CellIndex(:,j)==cellindex,J); % leaf nodes of this cell
        netPts       = find(ismember(gMRA.IniLabels,leafindex));   % points in the cell
        %%% indices of the points outside the cell
        netPtsRest   = setdiff((1:nAllPts),netPts);     % points outside the cell
        %% initialize the center
        center       = gMRA.Centers{cellindex};           % center of the cell
        rmax         = 0;
        TestPoint    = netPts;
        nTestPoint   = length(TestPoint);
        time         = Times;
        while time>0 && nTestPoint>0
            time = time-1;
            %%% distance between all points and the center
            DistAllPts   = sqrt(sum((bsxfun(@minus,X, center)).^2,1));            
            %%% radius test
            r  = gMRA.Radii(cellindex);
            %%% number of points in the cell within the r-ball
            netPtsIn   = find(DistAllPts(netPts)<=r);
            nPtsIn     = length(netPtsIn);
            nPtsInRest = length(find(DistAllPts(netPtsRest)<=r));
            ratio      = nPtsIn/nPtsInRest;
            while ratio < tol
                r = r/2;
                netPtsIn   = find(DistAllPts(netPts)<=r);
                nPtsIn     = length(netPtsIn);
                nPtsInRest = length(find(DistAllPts(netPtsRest)<=r));
                ratio      = nPtsIn/nPtsInRest;
            end
            rmax = max(rmax,r);
            %% reassign the center
            TestPoint   = setdiff(TestPoint,netPts(netPtsIn)); 
            nTestPoint  = length(TestPoint);  
            if nTestPoint>0
                centerindex = TestPoint(randi(nTestPoint));
                center      = X(:,centerindex);
            end
       end
        %%% record data
        BallInCell{j}(iter,1) = cellindex;
        BallInCell{j}(iter,2) = rmax;
        BallInCell{j}(iter,3) = length(netPts)/nAllPts;
    end
end

