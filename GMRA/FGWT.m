function XGWT = FGWT(gMRA, X)

%% FGWT - Forward Geometric Wavelet Transform on any input data
%
% INPUT
%   gMRA: the geometric multi-resolution analysis structure for X, as created by GMRA
%      X: D-by-N matrix of data points
%
% OUTPUT
%   XGWT: a structure of the following fields:
%      .leafNodeLabels  : an N-vector of indices of leaf nodes to which the input data are assigned by proximity
%      .LeafNodeSizes   : number of points assigned to each leaf node
%      .dists2NearestLeafNodes: corresponding distances from data points to the centers of the nearest leaf nodes
%      .CelWavCoeffs    : nLeafNodes by J matrix of cells, each cell (i,j) is a matrix containing in rows the wavelet coefficients of
%                           the points in the leaf node i and at the corresponding scale j
%      .CelScalCoeffs   : similar to above, but stores scaling coefficients
%      .CelTangCoeffs   : similar to above, but stores tangential coefficients
%      .MatWavCoeffs    : matrix of wavelet coefficients, rows correspond to scales, columns correspond to points
%      .maxWavDims      : vector of maximal wavelet dimension at each scale
%      .MatWavDims      : N-by-J matrix of wavelet dimensions
%      .CoeffsCosts     : cost in storing all nonzero wavelet coefficients

%% Parameters
J          = max(gMRA.Scales);                                                                                                       % Number of scales
nLeafNodes = numel(gMRA.LeafNodes);

%% Initialization
XGWT = struct();

% Subtract mean and project if needed
if isfield(gMRA.opts,'Mean') && ~isempty(gMRA.opts.Mean)
    X = bsxfun(@minus,X,gMRA.opts.Mean);
end
if isfield(gMRA.opts,'Proj') && ~isempty(gMRA.opts.Proj)
    X = gMRA.opts.Proj*X;
end

% Find the nearest leaf to each data point
if isfield( gMRA, 'TreeOnLeaves' )
    CoverTree = gMRA.TreeOnLeaves;
else
    CoverTree = [];
end

if gMRA.opts.verbose > 0, fprintf('\n\t FGWT: assigning points to leaves...'); end
XGWT.Timing.nrsearch                                = cputime;
[~, nearestLeafIdx, XGWT.dists2NearestLeafNodes]    = nrsearch(cat(2,gMRA.Centers{gMRA.LeafNodes}), X, 1, 0, struct('ReturnAsArrays',1,'NNInfo',struct('CoverTree',CoverTree)));
XGWT.Timing.nrsearch                                = cputime-XGWT.Timing.nrsearch;
if gMRA.opts.verbose > 0, fprintf('done. (%.5f secs)',XGWT.Timing.nrsearch); end

XGWT.leafNodeIdxs   = uint32(nearestLeafIdx);
XGWT.leafNodeSizes  = zeros(1,nLeafNodes);
XGWT.CelWavCoeffs   = cell(nLeafNodes,J);
XGWT.CelScalCoeffs  = cell(nLeafNodes,J);
XGWT.PointsInNet    = cell(length(gMRA.cp),1);
if gMRA.opts.addTangentialCorrections,
    XGWT.CelTangCoeffs = cell(nLeafNodes,J);
end
XGWT.Cel_cpidx      = zeros(nLeafNodes,J);
XGWT.DeltaPoint     = zeros(size(X,2),J);
XGWT.DeltaMatrix    = zeros(nLeafNodes,J);
XGWT.DeltaMatrixInf = zeros(nLeafNodes,J);
XGWT.CellRadii      = zeros(nLeafNodes,J);
XGWT.CellSizes      = zeros(nLeafNodes,J);

%% Pre-process for book-keeping and memory allocation and greater parallelism down the road
if gMRA.opts.verbose > 0, fprintf('\n\t FGWT: pre-processing...'); end
XGWT.Timing.preprocessing                           = cputime;

leafIdxswithPts                                     = unique(XGWT.leafNodeIdxs);
leafNodesWithPts                                    = gMRA.LeafNodes(leafIdxswithPts);
[sortedLeafNodeLabels,sortedLeafNodeLabels_idxs]    = sort(XGWT.leafNodeIdxs);
sortedLeafNodeLabels_idxs                           = uint32(sortedLeafNodeLabels_idxs);
sortedLeafNodeLabels_counts                         = uint32(histc(sortedLeafNodeLabels,unique(sortedLeafNodeLabels)));
sortedLeafNodeLabels_cumcounts                      = [0;uint32(cumsum(sortedLeafNodeLabels_counts))];
XGWT.leafNodeSizes                                  = sortedLeafNodeLabels_counts;

%% Count #points
PointsInNetCount                                    = zeros(1,length(gMRA.cp),'uint32');
PointsInNetCurIdx                                   = zeros(1,length(gMRA.cp),'uint32');
for i = 1:length(leafIdxswithPts)                                                                                               % Pre-processing for memory pre-allocation
    iFineNet                                        = leafNodesWithPts(i);
    PointsInNetCount(iFineNet)                      = PointsInNetCount(iFineNet) + XGWT.leafNodeSizes(i);
    if XGWT.leafNodeSizes(i)>0                                                                                                  % If not empty
        j = gMRA.Scales(iFineNet);                                                                                              % Current scale
        while j>1                                                                                                               % At the top of the tree
            iCoarseNet                              = gMRA.cp(iFineNet);                                                        % For the other scales, get the parent
            PointsInNetCount(iCoarseNet)            = PointsInNetCount(iCoarseNet) + XGWT.leafNodeSizes(i);
            iFineNet                                = iCoarseNet;
            j                                       = j-1;
        end
    end
end

%% Record points in the leaf nodes at all scales
for i = 1:length(leafIdxswithPts)                                                                                               % Pre-processing for memory allocation and indices
    iFineNet                                        = leafNodesWithPts(i);
    if PointsInNetCurIdx(iFineNet)==0                                                                                           % First time at this node:
        XGWT.PointsInNet{iFineNet}                  = zeros(1,PointsInNetCount(iFineNet),'uint32');                             %   allocate memory
        PointsInNetCurIdx(iFineNet)                 = 1;                                                                        %   mark node as visited
    end
    if XGWT.leafNodeSizes(i)>0                                                                                                  % If not empty
        idxs_tmp                                    = PointsInNetCurIdx(iFineNet):PointsInNetCurIdx(iFineNet)+XGWT.leafNodeSizes(i)-1;
        XGWT.PointsInNet{iFineNet}(idxs_tmp)        = sortedLeafNodeLabels_idxs(sortedLeafNodeLabels_cumcounts(i)+1:sortedLeafNodeLabels_cumcounts(i+1));
        PointsInNetCurIdx(iFineNet)                 = idxs_tmp(end)+1;
        j = gMRA.Scales(iFineNet);                                                                                              % Current scale
        while j>1                                                                                                               % At the top of the tree
            iCoarseNet = gMRA.cp(iFineNet);                                                                                     % For the other scales, get the parent
            if PointsInNetCurIdx(iCoarseNet)==0                                                                                 % First time at this node:
                XGWT.PointsInNet{iCoarseNet}        = zeros(1,PointsInNetCount(iCoarseNet),'uint32');                           %   allocate memory
                PointsInNetCurIdx(iCoarseNet)       = 1;                                                                        %   mark node as visited
            end
            idxs_tmp                                = PointsInNetCurIdx(iCoarseNet):PointsInNetCurIdx(iCoarseNet)+XGWT.leafNodeSizes(i)-1;
            XGWT.PointsInNet{iCoarseNet}(idxs_tmp)  = sortedLeafNodeLabels_idxs(sortedLeafNodeLabels_cumcounts(i)+1:sortedLeafNodeLabels_cumcounts(i+1));%XGWT.leafNodeSizes(i);
            PointsInNetCurIdx(iCoarseNet)           = idxs_tmp(end)+1;
            iFineNet                                = iCoarseNet;
            j                                       = j-1;
        end
    end
end
XGWT.Timing.preprocessing                            = cputime - XGWT.Timing.preprocessing;
if gMRA.opts.verbose > 0, fprintf('done. (%.5f secs)',XGWT.Timing.preprocessing); end


%% Process each leaf node separately and compute scaling and wavelet coefficients on the path to the root
if gMRA.opts.verbose > 0, fprintf('\n\t FGWT: computing coefficients...'); end
XGWT.Timing.coefficients                    = cputime;

leafIdxsnoPts                                       = setdiff(1:length(gMRA.LeafNodes),leafIdxswithPts);
if gMRA.opts.ComputeWavelets                                                                                                    %% Geometric wavelets are computed in GMRA
    for i = 1:length(leafIdxswithPts)                                                                                                
        iFineNet    = leafNodesWithPts(i);
        curIdx      = leafIdxswithPts(i);        
        pt_idxs                         = XGWT.PointsInNet{iFineNet};                                                           % Get the points in this leaf        
        j                               = gMRA.Scales(iFineNet);                                                                % Current scale
        if j==1                                                                                                                 % At the top of the tree
            XGWT.CelWavCoeffs{curIdx,1}  = gMRA.WavBases{iFineNet}*bsxfun(@minus, X(:,pt_idxs), gMRA.Centers{iFineNet});
            XGWT.Cel_cpidx(curIdx,1)     = iFineNet;
            XGWT.CelScalCoeffs{curIdx,1} = XGWT.CelWavCoeffs{curIdx,1};
            XGWT.CellRadii(curIdx,1)     = gMRA.Radii(iFineNet);
            XGWT.CellSizes(curIdx,1)     = PointsInNetCount(iFineNet);
        else
            iCoarseNet = gMRA.cp(iFineNet);                                                                                     % For the other scales, get the parent            
            if ~gMRA.opts.orthogonalizing && gMRA.opts.addTangentialCorrections
                XGWT.CelScalCoeffs{curIdx,j}    = gMRA.ScalBasis{iFineNet}*bsxfun(@minus, X(:,pt_idxs), gMRA.Centers{iFineNet});
                Projections_jmax                = X(:,pt_idxs);                                                                 % bsxfun(@plus, gMRA.ScalBasis{iFineNet}'*XGWT.CelScalCoeffs{curIdx,j}, gMRA.Centers{iFineNet});
            else
                Projections_jmax                = X(:,pt_idxs);
                XGWT.CelScalCoeffs{iFineNet,j}  = gMRA.ScalBasis{iFineNet}*bsxfun(@minus, Projections_jmax, gMRA.Centers{iFineNet});
            end                       
            if ~isempty(gMRA.WavBases{iFineNet})                
                    XGWT.CelWavCoeffs{curIdx,j}  = ComputeWaveletCoeffcients(    ...                                            % Compute wavelet coefficients
                                                                XGWT.CelScalCoeffs{curIdx,j}, gMRA.ScalBasis{iFineNet}, ...
                                                                gMRA.WavBases{iFineNet}, gMRA.opts.sparsifying, ...
                                                                gMRA.opts.precision, gMRA.O{iFineNet} );
            end            
            if gMRA.opts.addTangentialCorrections
                XGWT.CelTangCoeffs{curIdx,j}  = zeros(size(gMRA.ScalBasis{iCoarseNet},1),XGWT.leafNodeSizes(i));
            end
            
            XGWT.Cel_cpidx(curIdx,j)     = iFineNet;
            XGWT.CellRadii(curIdx,j)     = gMRA.Radii(iFineNet);   % radius
            XGWT.CellSizes(curIdx,j)     = PointsInNetCount(iFineNet);          
            % Coefficients in the wavelet direction: XGWT.CelWavCoeffs{i,j} + gMRA.WavConstsCoeffs{iFineNet} 
            % bsxfun(@plus,XGWT.CelWavCoeffs{i,j},gMRA.WavConstsCoeffs)
            % Tangential correction coefficients: XGWT.CelTangCoeffs{i,j}
            if ~isempty(gMRA.WavBases{iFineNet})  
                DeltaCoeffs  = [bsxfun(@plus,XGWT.CelWavCoeffs{curIdx,j},gMRA.WavConstsCoeffs{iFineNet} ) ; XGWT.CelTangCoeffs{curIdx,j}];
            else
                DeltaCoeffs  = XGWT.CelTangCoeffs{curIdx,j};
            end
            FineCoarseDiff2                 = sum(DeltaCoeffs.^2,1);
            XGWT.DeltaPoint(pt_idxs,j-1)    = FineCoarseDiff2;
            XGWT.DeltaMatrix(curIdx,j-1)    = sum(FineCoarseDiff2);
            XGWT.DeltaMatrixInf(curIdx,j-1) = max(sqrt(FineCoarseDiff2));
            
            j = j-1;                      
            
            while j>=1                                                                                                          % Go through scales between the leaf node and the root, from fine to coarse                                
                iFinerNet                   = iFineNet;
                iFineNet                    = iCoarseNet;
                iCoarseNet                  = gMRA.cp(iFineNet);
                XGWT.Cel_cpidx(curIdx,j)    = iFineNet;
                XGWT.CellRadii(curIdx,j)    = gMRA.Radii(iFineNet);   % radius
                XGWT.CellSizes(curIdx,j)    = PointsInNetCount(iFineNet);
                if ~isempty(gMRA.ScalBasisChange{iFinerNet}) && ~isempty(XGWT.CelScalCoeffs{curIdx,j+1})                   
                    if ~isempty(gMRA.ScalWavBasisChange{iFinerNet})                                                                                             %XGWT.CelScalCoeffs{i,j} =  bsxfun(@minus, gMRA.ScalBasis{iFineNet}*Projections, gMRA.ScalBasis{iFineNet}*gMRA.Centers{iFineNet});
                        Projections_compr       = bsxfun(@plus, gMRA.ScalBasisChange{iFinerNet}*XGWT.CelScalCoeffs{curIdx,j+1} - ...
                                                                gMRA.ScalWavBasisChange{iFinerNet}*XGWT.CelWavCoeffs{curIdx,j+1}, ...
                                                                gMRA.ScalCenterConsts{iFinerNet});                                                              
                    else                        
                        Projections_compr       = bsxfun(@plus, gMRA.ScalBasisChange{iFinerNet}*XGWT.CelScalCoeffs{curIdx,j+1}, ...
                                                                gMRA.ScalCenterConsts{iFinerNet});
                    end
                    if gMRA.opts.addTangentialCorrections && ~isempty(XGWT.CelTangCoeffs{curIdx,j+1})
                        Projections_compr       = Projections_compr - XGWT.CelTangCoeffs{curIdx,j+1};
                    end
                    XGWT.CelScalCoeffs{curIdx,j}     =  bsxfun(@minus,Projections_compr, gMRA.ScalCenter{iFineNet});
                elseif ~isempty(gMRA.ScalBasis{iFineNet})
                    XGWT.CelScalCoeffs{curIdx,j} = repmat(gMRA.ScalBasis{iFineNet}*(gMRA.Centers{iFinerNet}-gMRA.Centers{iFineNet}),1,length(XGWT.PointsInNet{iFinerNet}));
                end         
                
                if j==1, break;     end                                                                                        % Root is slightly different (wavelet subspace=scaling subspace)
                
                if ~(isempty(gMRA.WavBases{iFineNet}) || isempty(XGWT.CelScalCoeffs{curIdx,j})) % Wavelet coefficients
                        XGWT.CelWavCoeffs{curIdx,j}  = ComputeWaveletCoeffcients(    XGWT.CelScalCoeffs{curIdx,j}, gMRA.ScalBasis{iFineNet}, ...
                                                                                gMRA.WavBases{iFineNet}, gMRA.opts.sparsifying, ...
                                                                                gMRA.opts.threshold0 );
                end
                if gMRA.opts.addTangentialCorrections
                    if ~isempty(gMRA.ScalBasisChange{iFineNet}) && ~isempty(XGWT.CelScalCoeffs{curIdx,j})
                        XGWT.CelTangCoeffs{curIdx,j}     = bsxfun(@plus, gMRA.ScalBasisChange{iFineNet}*XGWT.CelScalCoeffs{curIdx,j}, ...     % Tangential corrections
                            gMRA.ScalBasisCoarserCenter{iFineNet}) - ...
                            gMRA.ScalBasis{iCoarseNet}*Projections_jmax;                                            %                    tangCoeffs                  = gMRA.ScalBasis{iCoarseNet}*(Projections-Projections_jmax);
                    end
                end    
            
                % Coefficients in the wavelet direction: XGWT.CelWavCoeffs{iLeaf,j} + gMRA.WavConstsCoeffs{iFineNet} 
                % bsxfun(@plus,XGWT.CelWavCoeffs{iLeaf,j},gMRA.WavConstsCoeffs)
                % Tangential correction coefficients: XGWT.CelTangCoeffs{iLeaf,j}
                if ~isempty(gMRA.WavBases{iFineNet}) && ~isempty(XGWT.CelTangCoeffs{curIdx,j})
                    DeltaCoeffs  = [bsxfun(@plus,XGWT.CelWavCoeffs{curIdx,j},gMRA.WavConstsCoeffs{iFineNet} ) ; XGWT.CelTangCoeffs{curIdx,j}];
                else
                    DeltaCoeffs  = XGWT.CelTangCoeffs{curIdx,j};
                end
                if ~isempty(DeltaCoeffs)
                    FineCoarseDiff2                = sum(DeltaCoeffs.^2,1);
                    XGWT.DeltaPoint(pt_idxs,j-1)   = FineCoarseDiff2;
                    XGWT.DeltaMatrix(curIdx,j-1)   = sum(FineCoarseDiff2);
                    XGWT.DeltaMatrixInf(curIdx,j-1)= max(sqrt(FineCoarseDiff2));
                end
                
                j = j-1;
            end
            
            XGWT.CelWavCoeffs{curIdx,1}      = XGWT.CelScalCoeffs{curIdx,1};
            if ~isempty(XGWT.CelWavCoeffs{curIdx,1})                                
                XGWT.CelTangCoeffs{curIdx,j} = bsxfun(@plus,XGWT.CelScalCoeffs{curIdx,j},gMRA.ScalCenter{iFineNet})-gMRA.ScalBasis{iFineNet}*Projections_jmax;                    %XGWT.CelTangCoeffs{i,j} = gMRA.ScalBasis{iFineNet}*(bsxfun(@plus,gMRA.ScalBasis{iFineNet}'*XGWT.CelScalCoeffs{i,j},gMRA.Centers{iFineNet})-Projections_jmax);               %XGWT.CelTangCoeffs{i,j} = gMRA.ScalBasis{iFineNet}*(Projections-Projections_jmax);
            end
        end
    end
else
    for i = 1:length(leafIdxswithPts)                                                                                                
        iFineNet    = leafNodesWithPts(i);
        curIdx      = leafIdxswithPts(i);
        %% leaf with points
        pt_idxs                         = XGWT.PointsInNet{iFineNet};                                                       % Get the points in this leaf        
        j                               = gMRA.Scales(iFineNet);                 
        % Scaling coeffients and projections at the current scale
        FineBasis         = gMRA.ScalBasis{iFineNet};
        FineScalCoeffs    = FineBasis*(bsxfun(@minus, X(:,pt_idxs), gMRA.Centers{iFineNet}));
        FineProjections   = bsxfun(@plus, FineBasis'*FineScalCoeffs, gMRA.Centers{iFineNet});
        Projections_jmax  = X(:,pt_idxs); %FineProjections;
        if j==1
            XGWT.CelScalCoeffs{curIdx,j} = FineScalCoeffs; 
            XGWT.Cel_cpidx(curIdx,j)     = iFineNet;  
            XGWT.CellRadii(curIdx,j)     = gMRA.Radii(iFineNet);
            XGWT.CellSizes(curIdx,j)     = PointsInNetCount(iFineNet);       
        else % j>1 record the current scale and pass it to the coarser scale
            while j>1
                % record the current scale
                XGWT.CelScalCoeffs{curIdx,j} = FineScalCoeffs; 
                XGWT.Cel_cpidx(curIdx,j)     = iFineNet;              
                XGWT.CellRadii(curIdx,j)     = gMRA.Radii(iFineNet);
                XGWT.CellSizes(curIdx,j)     = PointsInNetCount(iFineNet);               
                % look at the coarse scale
                iCoarseNet        = gMRA.cp(iFineNet); 
                CoarseBasis       = gMRA.ScalBasis{iCoarseNet};
                CoarseScalCoeffs  = CoarseBasis*(bsxfun(@minus, Projections_jmax, gMRA.Centers{iCoarseNet}));
                CoarseProjections = bsxfun(@plus, CoarseBasis'*CoarseScalCoeffs, gMRA.Centers{iCoarseNet});
                
                % compute Delta
                FineCoarseDiff2 = sum((FineProjections - CoarseProjections).^2,1);
                
                XGWT.DeltaPoint(pt_idxs,j-1)      = FineCoarseDiff2;
                XGWT.DeltaMatrix(curIdx,j-1)      = sum(FineCoarseDiff2);
                XGWT.DeltaMatrixInf(curIdx,j-1)   = max(sqrt(FineCoarseDiff2));
                % pass to the next scale
                j = j-1;
            
                iFineNet   = iCoarseNet;
            
                FineScalCoeffs  = CoarseScalCoeffs;  %bsxfun(@minus, X(:,pt_idxs), gMRA.Centers{iFineNet})'*FineBasis;
                FineProjections = CoarseProjections; %bsxfun(@plus, FineBasis*FineScalCoeffs', gMRA.Centers{iFineNet});
          
            end
            % j = 1
            FineBasis       = gMRA.ScalBasis{iFineNet};
            FineScalCoeffs  = FineBasis*(bsxfun(@minus, Projections_jmax, gMRA.Centers{iFineNet}));
        
            XGWT.CelScalCoeffs{curIdx,j} = FineScalCoeffs; 
            XGWT.Cel_cpidx(curIdx,j)     = iFineNet;         
            XGWT.CellRadii(curIdx,j)     = gMRA.Radii(iFineNet);
            XGWT.CellSizes(curIdx,j)     = PointsInNetCount(iFineNet);           
        end
   end
end

%% leaf without points
for i = 1:length(leafIdxsnoPts)
    iFineNet    = gMRA.LeafNodes(leafIdxsnoPts(i));
    curIdx      = leafIdxsnoPts(i);
    j           = gMRA.Scales(iFineNet);                                                                                              % Current scale
    if j==1                                                                                                                 % Single-node tree
        XGWT.Cel_cpidx(curIdx,1)     = iFineNet;
        XGWT.CellRadii(curIdx,1)     = gMRA.Radii(iFineNet);
        XGWT.CellSizes(curIdx,1)     = PointsInNetCount(iFineNet);           
    else
        iCoarseNet               = gMRA.cp(iFineNet);                                                                                     % For scale at current leaf node
        XGWT.Cel_cpidx(curIdx,j) = iFineNet;
        j = j-1;
        while j>1                                                                                                           % For scales between the leaf node and the root
            iFineNet = iCoarseNet;
            iCoarseNet = gMRA.cp(iFineNet);
            XGWT.Cel_cpidx(curIdx,j)     = iFineNet;
            XGWT.CellRadii(curIdx,j)     = gMRA.Radii(iFineNet);   % radius
            XGWT.CellSizes(curIdx,j)     = PointsInNetCount(iFineNet);           
            j = j-1;
        end
            
        XGWT.Cel_cpidx(curIdx,1)     = iCoarseNet;
        XGWT.CellRadii(curIdx,1)     = gMRA.Radii(iCoarseNet);   % radius
        XGWT.CellSizes(curIdx,1)     = PointsInNetCount(iCoarseNet);            
    end
    clearvars XGWT.CelWavCoeffs XGWT.CelTangCoeffs
end

XGWT.Timing.coefficients                             = cputime-XGWT.Timing.coefficients;
if gMRA.opts.verbose > 0, fprintf('done. (%.5f secs)',XGWT.Timing.coefficients); end

%% Compute cell-wise Deltajk
XGWT.DeltaCell    = zeros(nLeafNodes,J); 
XGWT.DeltaCellInf = zeros(nLeafNodes,J);    
for j = 1 : J
    UniqueIndexj    = unique(XGWT.Cel_cpidx(:,j));
    UniqueIndexj    = UniqueIndexj(UniqueIndexj>0);
    nUniqueIndexj   = numel(UniqueIndexj);
    XGWTCel_cpidxj  = XGWT.Cel_cpidx(:,j);
    for i = 1 : nUniqueIndexj
        currentindex                  = UniqueIndexj(i);
        rowindex                      = find(XGWTCel_cpidxj==currentindex);
        currentdelta                  = sum (XGWT.DeltaMatrix(rowindex,j))/size(X,2);%/currentnPts;
        XGWT.DeltaCell(rowindex,j)    = sqrt(currentdelta);
        XGWT.DeltaCellInf(rowindex,j) = max(XGWT.DeltaMatrix(rowindex,j));
    end
end
XGWT.DeltaPoint  = sqrt(XGWT.DeltaPoint);
XGWT.CellMeasure = XGWT.CellSizes/size(X,2);

%% Compute coefficient cost
XGWT.CoeffsCosts = sum(sum(cellfun(@(x) numel(x),XGWT.CelWavCoeffs)));

return


function wavCoeffs = ComputeWaveletCoeffcients(data_coeffs, scalBases, wavBases, sparsifying, epsilon, transform)

if ~sparsifying,
    if nargin<6 || isempty(transform)
        wavCoeffs = wavBases*scalBases'*data_coeffs;
    else
        data_transformed = bsxfun(@plu,transform.T'*transform.b * scalBases'*data_coeffs,transform.c');
        wavCoeffs = wavBases*data_transformed;
    end
else
    %wavCoeffs = transpose(mexLasso(data,wavBases,struct('iter', -5, 'mode', 1, 'lambda', 1e-3)));
    %wavCoeffs = full(omp(wavBases'*data,wavBases'*wavBases,5))';
    wavCoeffs = full(omp2(data*wavBases,sum(data.*data),wavBases*wavBases',epsilon))';
end

return;
