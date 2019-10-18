function DataAGWT = FAGWT_TruncateTree(XGWT,kappa)
%% Perform adaptive GMRA
%% Step 1: Threshold the tree for a specific tau
%%         Refinement criterion is based on training samples
%%  XGWT:  It contains the refinement criterion: Deltajk
%%  kappa: Thresholding parameter in adaptive GMRA

%% threshold
DataAGWT.kappa        = kappa;
n                     = size(XGWT.DeltaPoint,1);
nLeafNodes            = size(XGWT.Cel_cpidx,1);
DataAGWT.tau          = kappa*sqrt(log(n)/n);
DataAGWT.ThreshKappa  = XGWT.CellRadii*DataAGWT.tau;


%% Truncate the tree and find adaptive paritition
DataAGWT.DT             = XGWT.DeltaCell > DataAGWT.ThreshKappa;
DataAGWT.PartitionScale = zeros(nLeafNodes,1);
DataAGWT.PartitionIndex = zeros(nLeafNodes,1);

for i = 1 : nLeafNodes
    LastOne = find(DataAGWT.DT(i,:),1,'last');
    if isempty(LastOne)
        LastOne = 0;
    end
    LastOne = LastOne + 1; % take the outer leaves
    while LastOne>0 && XGWT.Cel_cpidx(i,LastOne) == 0
        LastOne = LastOne - 1;
    end
    if LastOne == 0,    LastOne = 1;    end;
    DataAGWT.PartitionScale(i) = LastOne;
    DataAGWT.PartitionIndex(i) = XGWT.Cel_cpidx(i,LastOne);
end
    
return
