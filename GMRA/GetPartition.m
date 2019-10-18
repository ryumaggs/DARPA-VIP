function Partition = GetPartition( gMRA, XGWT, opts )

%
% function Partition = GetPartition( gMRA, XGWT, opts )
%
% IN:
%   gMRA        : GMRA data structure
%   XGWT        : FGWT of data
%   [opts]      : data structure of options:
%       [type]  : value in {'Uniform','Adaptive','OuterLeaves'}. Default: 'Uniform'.
%       [scale] : scale of partition if type=='Uniform'. Defaults to reasonable value.
%       [kappa] : threshold determining the partition if type=='Adaptive'. Defaults to reasonable value.
%       [MarkedNodes] : constructs partition out of outer leaves of the smallest tree containing the nodes MarkedNodes, if type=='OuterLeaves'.
%       [disregsardleaves] : disregards the leaves of the GMRA tree when computing the adaptive partition. Default: false.
%       [sigma] : noise level
%       [refinementnorm]: norm of refinement criterion
%       [ifscaled]: 1: scale-dependent threshold 0: scale-independent
%
% OUT:
%   Partition   : structure with the following fields:
%       leafParentInPartition : vector of length equal to the number of leafnodes in the GMRA tree, with k-th entry being
%                               the index of the node in the GMRA tree to which that leafnode belongs in the partition
%       nodes                 : list of nodes in the GMRA tree associated with the partition. Note that these are not disjoint,
%                               they are simply unique(leafParentInPartition)
%
% (c) Mauro Maggioni
%

if nargin<3,     opts = [];                             end
if ~isfield(opts,'type'),   opts.type = 'Uniform';      end
if ~isfield(opts,'sigma'),  opts.sigma= 0;      end
if ~isfield(opts,'refinementnorm'), opts.refinementnorm = '2'; end
if ~isfield(opts,'ifscaled'), opts.ifscaled = 1; end

nLeafNodes   = length(gMRA.LeafNodes);

if strcmpi(opts.type,'Uniform')                                                                                                 %% Uniform partition
    if ~isfield(opts,'scale') || isempty(opts.scale)
        opts.scale = floor(mean(gMRA.Scales));
    end
    opts.scale = min(opts.scale,find(gMRA.Radii_AvePerScale>=opts.sigma,1,'last')); % stop above the level of noise
    %Partition = sort([find(gMRA.Scales == opts.scale);gMRA.LeafNodes(gMRA.Scales(gMRA.LeafNodes)<opts.scale)], 'ascend');      % DEBUG only
    Partition.leafParentInPartition = zeros(nLeafNodes,1);
    Partition.leafScale             = zeros(nLeafNodes,1);
    for i = 1 : nLeafNodes                                                                                                      % For every leaf node, find nontrivial node at finest scale, but no larger than j
        LastOne = opts.scale; %min(...,find(XGWT.CellSizes(i,:)>=opts.minleafsize,1,'last'));
        while LastOne>0 && XGWT.Cel_cpidx(i,LastOne)==0
            LastOne = LastOne - 1;
        end
        if LastOne == 0,        LastOne = 1;    end
        Partition.leafParentInPartition(i) = max(XGWT.Cel_cpidx(i,LastOne),1);                                                  % Save index of selected node, for every leaf node
        Partition.leafScale(i)             = LastOne;
    end
    Partition.nodes = unique(Partition.leafParentInPartition);
elseif strcmpi(opts.type,'Adaptive')                                                                                            %% Adaptive partition
    if ~isfield(opts,'kappa')
        opts.kappa = [];
    end
    if ~isfield(opts,'disregardleaves')
        opts.disregardleaves = false;
    end
    
    n   = size(XGWT.DeltaPoint,1);                                                                                              % Set threshold
    if ~isempty(opts.kappa)
        tau          = opts.kappa*sqrt(log(n)/n);
    else
        tau          = median(abs(XGWT.DeltaCell))/median(XGWT.CellRadii);
    end
    SubtreeIndicator                = XGWT.CellRadii>opts.sigma; %(XGWT.CellSizes>=opts.minleafsize).*(XGWT.CellRadii>opts.sigma); % partition is taken on the subtree above the noise level whose leaf size >= opts.minleafsize
    for i = 1:nLeafNodes
        SubtreeIndicator(i,find(SubtreeIndicator(i,:),1,'last')) = 0;      
    end
    % threshold
    if opts.ifscaled
        ThreshKappa                     = XGWT.CellRadii*tau;
    else
        ThreshKappa                     = tau;
    end
    % threshold the tree 
   switch opts.refinementnorm
        case '2'
            DeltaCell                       = XGWT.DeltaCell;
        case 'inf'
            DeltaCell                       = XGWT.DeltaCellInf;
   end
   if opts.disregardleaves
        DeltaCell(:,end-2:end)      = 0;
    end
    DT                              = (DeltaCell> ThreshKappa).*SubtreeIndicator ;                                                             % Find wavelet energies above threshold
    Partition.leafParentInPartition = zeros(nLeafNodes,1);
    Partition.leafScale             = zeros(nLeafNodes,1);
    for i = 1 : nLeafNodes                                                                                                      % For every leaf node, find finest scale with wavelet energy above threshold
        LastOne = find(DT(i,:),1,'last');
        if isempty(LastOne),    LastOne = 1;    end
        LastOne = LastOne + 1; %take ourter leaves                                                                                                 % Take the outer leaves
%         if opts.disregardleaves
%             if XGWT.Cel_cpidx(i,LastOne)>0
%                 if gMRA.isaleaf( XGWT.Cel_cpidx(i,LastOne) )
%                     if LastOne > 1,
%                         LastOne = LastOne - 1;
%                     end
%                 end
%             end
%         end
        while LastOne>0 && XGWT.Cel_cpidx(i,LastOne) == 0 
            LastOne = LastOne - 1;
        end
        if LastOne == 0,        LastOne = 1;    end
        % make sure to choose a partition
        while ~isempty(find(DT(XGWT.Cel_cpidx(:,LastOne)==XGWT.Cel_cpidx(i,LastOne),LastOne:end),1)) 
            LastOne = LastOne+1;
        end 
        %
        Partition.leafParentInPartition(i) = max(XGWT.Cel_cpidx(i,LastOne),1);                                                  % Save index of outer leaves of adaptive tree, one for each leafnode
        Partition.leafScale(i)             = LastOne;
    end
    Partition.nodes = unique(Partition.leafParentInPartition);
elseif strcmpi(opts.type,'OuterLeaves')                                                                                            %% Adaptive partition
    DT = zeros(size(XGWT.Cel_cpidx));
    for k = 1:length(opts.MarkedNodes)
        [idxs_i,idxs_j] = find(XGWT.Cel_cpidx==tableMarkedNodes(k));
        for l = 1:length(idxs_i)
            DT(idxs_i(l),idxs_j(l)) = 1;
        end
    end
    Partition.leafParentInPartition = zeros(nLeafNodes,1);
    Partition.leafScale             = zeros(nLeafNodes,1);
    for i = 1 : nLeafNodes                                                                                                      % For every leaf node, find finest scale with wavelet energy above threshold
        LastOne = find(DT(i,:),1,'last');
        if isempty(LastOne),    LastOne = 1;    end
        LastOne = LastOne + 1; %take ourter leaves                                                                                                 % Take the outer leaves
        if opts.disregardleaves
            if XGWT.Cel_cpidx(i,LastOne)>0
                if gMRA.isaleaf( XGWT.Cel_cpidx(i,LastOne) )
                    if LastOne > 1,
                        LastOne = LastOne - 1;
                    end
                end
            end
        end
        while LastOne>0 && XGWT.Cel_cpidx(i,LastOne) == 0
            LastOne = LastOne - 1;
        end
        if LastOne == 0,        LastOne = 1;    end
        % make sure to choose a partition
        while ~isempty(find(DT(XGWT.Cel_cpidx(:,LastOne)==XGWT.Cel_cpidx(i,LastOne),LastOne:end),1)) 
            LastOne = LastOne+1;
        end 
        %
        Partition.leafParentInPartition(i) = max(XGWT.Cel_cpidx(i,LastOne),1);                                                  % Save index of outer leaves of adaptive tree, one for each leafnode
        Partition.leafScale(i)             = LastOne;
    end
    Partition.nodes = unique(Partition.leafParentInPartition);
end

return