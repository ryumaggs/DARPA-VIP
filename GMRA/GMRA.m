function gMRA = GMRA(X,opts)
%
% Geometric MultiResolution Analysis (GMRA) for Data Sets
%
% Usage: gMRA = GMRA(X, opts)
%
% INPUT:
%   X: D-by-N data matrix with each column being a data point
%   opts: a structure with the following optional parameters:
%       [precision]        : precision parameter for the leaf nodes; if provided, it will overwrite threshold0 at the leaf nodes. Default: 1e-2.
%       [ManifoldDimension]: if a positive number, then it will be provided to all nodes as local dimension,
%                            except for the leaf nodes whose dimensions are determined to achieve the requested precision.
%                            If zero (default), then the local dimensions at all nodes are determined adapatively. Default: 0.
%       [errorType]        : 'relative' or 'absolute' error computed from the local singular values. Default: 'relative'.
%       [threshold0]       : corresponding threshold to each error type, can be a scalar or a vector indexed from coarsest scale (1) to finest scale (J).
%                            (if length < J, the last element will be replicated for scales (J-length+1):J).
%                            Default: 0.5 (uniformly for all J scales). Used only if ManifoldDimension==0
%       [threshold1]       : threshold for determining the wavelet dimensions. Default: 1e-5.
%       [GWTversion]       : one of the following:
%                               = 0, plain construction
%                               = 1, orthogonal geometric wavelets
%                               = 2, minimal encoding-cost pruning
%                            Default: 0.
%       [addTangentialCorrections]: whether to add tangential corrections so as to use the best approximations at each scale. Default: true.
%       [PartitionType]    : a string indicating the type of partioning method to be used:
%                               'nesdis'    : dyadic nested dissection (metis). Default.
%                               'covertree' : overlapping balls at the appropriate scale, using covertrees.
%           [smallestMetisNet]    : lower bound for the sizes of the METIS nets. Default: 10.
%           [CoverTreeOpts]       : various options for covertree:
%                                   [MaxSamples]        : maximum number of points to feed to covertree. Default: Inf
%                                   [RefineCoverTree]   : refines the covertree by reshuffling points by the level in the non-refined covertree. Default: false
%           [CoverTreeBuildOpts]  : options for covertree construction. Passed to covertree_build (see there for details).
%           [CoverTreeExpandOpts] : options for expanding the covertree. Pass to covertree_expand (see there for details).
%           [CoverTreeTrimOpts]   : options for trimming the covertree. Pass to covertree_trim (see there for details).
%       [ComputeWavelets]  : Computes wavelets or not. It can be set to false only if GWTversion==0. Default: true.
%       [CoverTree]        : covertree. Uses this if provided.
%       [Predictor]        : structure with the following options:
%                               [PredictorType] : one of {'none','orthogonal'}. Default: none. MM:TBD
%
%   GRAPH OPTIONS - these options affect the construction of the graph, when a graph is constructed (either upon request or because METIS partitioning is requested)
%       [ConstructGraph]   : constructs graph (output in gMRA.Graph) or not. Default: 0, unless PartitionType=='nesdis'.
%       [knn]              : size of the neighborhood graph. Default: 50.
%       [knnAutotune]      : local scaling parameter. Default: 20.
%
%   TREE OPTIMIZATION OPTIONS
%       [sparsifying]      : whether to sparsify the wavelet bases using algorithms such as k-svd and SPAMS.
%                            Default: false.
%           [sparsifying_method]:'ksvd' (by M. Elad et al.) or 'spams' (by J. Mairal). Default: 'ksvd'.
%           [sparsifying_oversampling] : when sparsifying, will seek a dictionary of size (sparsifying_oversampling)*(dimension of wavelet subspace). Default: 2.
%       [splitting]        : whether to separate out the common intersection of the wavelet subspaces
%                            associated to the children of a fixed node in order to save the encoding cost.
%                            Default: false.
%       [threshold2]       : threshold for determining the intersection dimension of the wavelet subspaces associated to all children of a fixed node.
%                            Default: 0.1.
%       [mergePsiCapIntoPhi]:whether to merge the common part of the wavelet subspaces associated to the children into the scaling function of the parent.
%                            Default: false.
%       [shrinkage]        : 'soft' or 'hard' shrinkage of the wavelet coefficients. Default: 'hard'.
%       [coeffs_threshold] : threshold for shrinking wavelet coefficients. Default: 1e-10.
%       [verbose]          : level of verbosity, a number in {0,1}. Default: 0.
%       [MaxDim]           : maximum subspace dimension to be ever considered. Default:50.
%       [graph]            : graph structure as returned by GraphDiffusion. If provided, does not recompute the graph.
%       [tree]             : tree information as returned by nesdis (nested dissection). It will be gMRA.cp as described below. If provided, does not rerun nesdis or covertrees.
%
%   RUNNING OPTIONS
%       [parallel]         : run in parallel. Default: false.
%
%
% OUTPUT:
%   gMRA: a structure with the following fields:
%       .cp                : vector encoding the metis tree structure with .cp(x) = the parent of x. Each node is a subset of X and a parent is the union of its children
%       .LeafNodes         : vector of labels of the leaf nodes.
%       .isaleaf           : vector of 1 (leaf node) or 0 (non-leaf node).
%       .IniLabels         : the labels of the data points w.r.t. the leaf nodes in the tree (all the leaves are a full partition of the data).
%       .PointsInNet       : cell array, PointsInNet{i} contains the labels of all the points in the node i of the metis tree.
%       .Sizes             : vector of number of points in each node.
%       .Radii             : vector of the radius of each node.
%       .Scales            : vector of scales of all the nodes in the metis tree; the root has scale 1 and a larger scale implies a finer approximation.
%       .Centers           : cell array of the centers of the nodes.
%       .ScalBasis         : cell array of the local bases of the nodes.
%       .Sigmas            : cell array of the local singular values
%       .WavBasis          : cell array of wavelet bases; each cell encodes the basis vectors that are present in the current node but orthogonal to the
%                               parent. For simplicity, the wavelet basis at the root is identified with the scaling function at the root.
%       .WavConsts         : cell array of wavelet translations that are needed to move to a node from its parent.
%       .WavSingVals       : cell array of corresponding singular values associated with the wavelet bases. At the root, it coincides with .Sigmas
%       .epsEncodingCosts  : vector of total encoding costs at each node when approximating the local data by the lowest-dimensional pca plane within the given precision .threshold0.
%       .DictCosts         : overall cost for storing the dictionary (i.e., wavelet bases and constants). When .addTangentialCorrections = true, the scaling functions will also be included in the dictionary.
%       .Timing            : structure with the following fields:
%                               graph  : cputime (in seconds) taken by the graph construction
%                               nesdis : cputime (in seconds) taken by multiscale partitioning
%                               GW     : cputime (in seconds) for the whole construction
%
% Required Packages:
%   1. Diffusion Geometry [by Mauro Maggioni]
%   2. Metis              [by George Karypis et al.]
%   3. SuiteSparse        [by Tim Davis, for the metis and nesdis wrappers]
%   4. LightSpeed         [by Tom Minka]
% If sparsifying wavelet basis, then also need the following two packages:
%   5. K_SVD              [by Michael Elad et al.] and/or
%   6. SPAMS              [by Julien Mairal et al.]
%
% Publications:
%   1. Multiscale Geometric Methods for Data Sets II: Geometric Wavelets, W.K. Allard, G. Chen and M. Maggioni, ACHA, 2011
%   2. Multiscale Geometric Dictionaries for Point-Cloud Data, G. Chen, and M. Maggioni, The 9th International Conference on Sampling Theory and Applications (SampTA), Singapore, 2011
%   3. Multiscale Geometric Wavelets for the Analysis of Point Clouds, G. Chen and M. Maggioni, The 44th Annual Conference on Information Sciences and Systems (CISS), Princeton, NJ, 2010
%
% (c) 2011, Mauro Maggioni, Duke University
% Contact: mauro}@math.duke.edu

%clear functions

% %% Handle parallelism carefully
% % Close all existing Matlab parallel pools
% fprintf('\n');
% delete(gcp('nocreate'));
% fprintf('\n');

rng('default');
rng('shuffle');

%% Parameters
gMRA.Timing.GW = cputime;

if nargin<1,                            error('The input X must be provided.');                                                 end
if nargin<2,                            opts = struct();                                                                        end
if ~isfield(opts, 'GWTversion') || ...
        isempty(opts.GWTversion),        opts.GWTversion = 0;                                                                   end
opts.orthogonalizing  = (opts.GWTversion == 1);
opts.pruning          = (opts.GWTversion == 2);
if ~isfield(opts, 'ComputeWavelets'),    opts.ComputeWavelets = true;                                                           end

% Parameters for building neighborhood graph
if ~isfield(opts, 'smallestMetisNet'),   opts.smallestMetisNet = 30;                                                            end
if ~isfield(opts, 'PartitionType'),      opts.PartitionType = 'nesdis';                                                         end
if ~isfield(opts, 'CoverTreeOpts'),      opts.CoverTreeOpts = [];                                                               end
if ~isfield(opts.CoverTreeOpts,'MaxSamples'), opts.CoverTreeOpts.MaxSamples = Inf;                                              end
if ~isfield(opts.CoverTreeOpts,'RefineCoverTree'), opts.CoverTreeOpts.RefineCoverTree = false;                                  end
if ~isfield(opts, 'CoverTreeBuildOpts'), opts.CoverTreeBuildOpts = struct('theta',0.5,'numlevels',int32(1000),'minlevel',int32(0));  end
if ~isfield(opts, 'CoverTreeExpandOpts'),opts.CoverTreeExpandOpts = struct('ExtRange','max');                                   end
if ~isfield(opts, 'CoverTreeTrimOpts'),  opts.CoverTreeTrimOpts = struct('TrimType','Size','TrimValue',50 );                    end
if ~isfield(opts, 'ConstructGraph')      opts.ConstructGraph = false;                                                           end
if ~isfield(opts, 'knn'),                opts.knn = 50;                                                                         end
if ~isfield(opts, 'knnAutotune'),        opts.knnAutotune = 30;                                                                 end
if ~isfield(opts, 'parallel'),           opts.parallel = false;                                                                 end
if ~isfield(opts, 'Predictor'),          opts.Predictor = struct();                                                             end

% Parameters for choosing local PCA dimensions
if ~isfield(opts, 'ManifoldDimension'),  opts.ManifoldDimension = 0;
elseif opts.pruning && opts.ManifoldDimension>0, % conflict
    opts.ManifoldDimension = 0;
    warning('Manifold Dimension is NOT used and has been set to zero in order to allow for locally adaptive dimensions!'); %#ok<WNTAG>
end
if ~isfield(opts, 'threshold0'),                        opts.threshold0 = 0.5;                                                  end
if ~isfield(opts, 'errorType'),                         opts.errorType = 'relative';                                            end
if ~isfield(opts, 'precision'),                         opts.precision = 1e-2;                                                  end

switch opts.GWTversion
    case 0
        if ~isfield(opts, 'addTangentialCorrections'),  opts.addTangentialCorrections   = true;                                 end
        if ~isfield(opts, 'splitting'),                 opts.splitting                  = false;                                end
        if ~isfield(opts, 'sparsifying'),               opts.sparsifying                = false;                                end
        if ~isfield(opts, 'sparsifying_oversampling'),  opts.sparsifying_oversampling   = 2;                                    end
        if isfield(opts,  'avoidLeafnodePhi') warning(sprintf('\n\t GMRA:opts.avoidLeafnodePhi obsolete.'));                    end
        if opts.sparsifying && ...
                ~isfield(opts, 'sparsifying_method')    opts.sparsifying_method         = 'ksvd';                               end
    case {1,2}
        opts.addTangentialCorrections   = false;
        opts.splitting                  = false;
        if (opts.GWTversion == 2) || (~isfield(opts, 'sparsifying')),    opts.sparsifying = false;                              end
end
if ~isfield(opts, 'threshold1'),                        opts.threshold1                 = 1e-5;                                 end
if (opts.splitting || opts.pruning) && ...
        ~isfield(opts, 'threshold2'),                   opts.threshold2                 = 1e-1;                                 end
if opts.splitting && ...
        ~isfield(opts, 'mergePsiCapIntoPhi'),           opts.mergePsiCapIntoPhi = false;                                        end
if ~isfield(opts, 'shrinkage'),                         opts.shrinkage          = 'hard';                                       end
if ~isfield(opts, 'coeffs_threshold'),                  opts.coeffs_threshold   = 1e-10;                                        end
if ~isfield(opts, 'verbose'),                           opts.verbose            = 0;                                            end
if ~isfield(opts, 'MaxDim'),                            opts.MaxDim             = 50;                                           end
if ~isfield(opts, 'graph'),                             opts.graph              = [];                                           end
if ~isfield(opts, 'graphNormalization'),                opts.graphNormalization = 'beltrami';                                   end
if isfield(opts,'Proj') && ~isempty(opts.Proj)
    opts.AmbientDimension  = size(opts.Proj,2);
else
    opts.AmbientDimension  = size(X,1);
end
opts.EmbedDim = size(X,1);
opts.N        = size(X,2);

if ~isfield(opts.Predictor,'PredictorType'),    opts.Predictor.PredictorType = 'none';                                          end

%%
gMRA.opts   = opts;

%% Tree construction
if isempty(opts.graph),
    GraphDiffusionOpts          = struct('KNN', opts.knn,'kNNAutotune',opts.knnAutotune,'Normalization',opts.graphNormalization, 'Display',0,'kEigenVecs',50,'Symmetrization','W+Wt','NNMaxDim',0);
    if strcmpi( opts.PartitionType,'covertree' ),
        GraphDiffusionOpts.FastNNSearcher = 'covertree';
        % Construct the cover tree
        if ~isfield(opts,'CoverTree') || isempty(opts.CoverTree)
            if gMRA.opts.verbose>0; fprintf('\n Constructing covertree...'); end
            if opts.CoverTreeOpts.MaxSamples>=opts.N
                opts.CoverTreeOpts.SamplesIdxs  = [];
                opts.CoverTreeOpts.Samples      = X;                
            else
                opts.CoverTreeOpts.SamplesIdxs  = randperm(opts.N,opts.CoverTreeOpts.MaxSamples);
                opts.CoverTreeOpts.Samples      = X(:,opts.CoverTreeOpts.SamplesIdxs);
            end
            gMRA.Timing.covertree           = cputime;
            gMRA.CoverTree                  = covertree_build( opts.CoverTreeOpts.Samples,opts.CoverTreeBuildOpts );
            gMRA.CoverTree
            gMRA.Timing.covertree           = cputime-gMRA.Timing.covertree;
            if gMRA.opts.verbose>0; fprintf('done. (%.5f secs)',gMRA.Timing.covertree); end
        else
            gMRA.CoverTree = opts.CoverTree;
            gMRA.Timing.covertree = 0;
        end
        GraphDiffusionOpts.DistInfo.CoverTree = gMRA.CoverTree;
        if opts.CoverTreeOpts.RefineCoverTree
            if gMRA.opts.verbose>0; fprintf('\n Constructing refined covertree...'); end
            gMRA.Timing.covertree2  = cputime;
            gMRA.sorted_idxs        = covertree_sortPtsbyLevel( gMRA.CoverTree );
            if isempty(opts.CoverTreeOpts.Samples)
                gMRA.CoverTree      = covertree_build( X(:,gMRA.sorted_idxs),opts.CoverTreeBuildOpts );
            else
                gMRA.CoverTree      = covertree_build( opts.CoverTreeOpts.Samples(:,gMRA.sorted_idxs),opts.CoverTreeBuildOpts );
            end
            gMRA.Timing.covertree2  = cputime-gMRA.Timing.covertree2;
            if gMRA.opts.verbose>0; fprintf('done. (%.5f secs)',gMRA.Timing.covertree); end
        else
            gMRA.sorted_idxs        = [];
        end
    else
        GraphDiffusionOpts.FastNNSearcher = '';
    end
    if opts.ConstructGraph || strcmp(opts.PartitionType,'nesdis')
        if gMRA.opts.verbose,     fprintf('\n Constructing graph...'); end
        gMRA.Timing.graph       = cputime;
        try
            gMRA.Graph          = GraphDiffusion( X, 0, GraphDiffusionOpts );
        catch
            keyboard
        end
        gMRA.Timing.graph       = cputime-gMRA.Timing.graph;
        if gMRA.opts.verbose,     fprintf('done. (%.3f sec)',gMRA.Timing.graph); end
    end
end

%% Build the tree of multiscale partitions
if gMRA.opts.verbose,   fprintf('\n Constructing multiscale partitions of type %s...',opts.PartitionType); end; tic
switch opts.PartitionType,
    case 'nesdis',
        gMRA.Timing.nesdis          = cputime;
        [~,gMRA.cp,cmember]         = nesdis(gMRA.Graph.W,'sym',opts.smallestMetisNet);
        gMRA.isaleaf                = leafnodes(gMRA.cp)';
        gMRA.LeafNodes              = find(gMRA.isaleaf);
        gMRA.IniLabels              = dissolve_separator(X, gMRA.cp, cmember, gMRA.LeafNodes, 'centers');
        gMRA.Timing.nesdis          = cputime-gMRA.Timing.nesdis;
        gMRA.pc                     = cp2pc( gMRA.cp );
    case 'covertree'
        if gMRA.opts.verbose>0; fprintf('\n\t First trim to the covertree...'); end
        gMRA.Timing.covertree_firsttrim = cputime;
        [~,cp_idxs,~,TreeTrimmed]       = covertree_trim( gMRA.CoverTree, opts.CoverTreeTrimOpts );                              % Trim the tree
        gMRA.Timing.covertree_firsttrim = cputime-gMRA.Timing.covertree_firsttrim;
        if gMRA.opts.verbose>0; fprintf('done. (%.3f sec)',gMRA.Timing.covertree_firsttrim); end
        if gMRA.opts.verbose>0; fprintf('\n\t Expanding the covertree...'); end
        gMRA.Timing.covertree_expand    = cputime;
        [ExpTree,gMRA.cp_nIdxs]         = covertree_expand( TreeTrimmed, opts.CoverTreeExpandOpts );                            % Expand the covertree to a more regular tree
        gMRA.Timing.covertree_expand    = cputime-gMRA.Timing.covertree_expand;
        if gMRA.opts.verbose>0; fprintf('done. (%.3f sec)',gMRA.Timing.covertree_expand); end
        clear TreeTrimmed
        if gMRA.opts.verbose>0; fprintf('\n\t Second trim to the covertree...'); end
        gMRA.Timing.covertree_trim      = cputime;
        [gMRA.cp,idxs_kept,gMRA.pc]     = covertree_trim( ExpTree, opts.CoverTreeTrimOpts );                                    % Trim back the tree to desired level
        gMRA.Timing.covertree_trim      = cputime-gMRA.Timing.covertree_trim;
        if gMRA.opts.verbose>0; fprintf('done. (%.3f sec)',gMRA.Timing.covertree_trim); end
        gMRA.cp_nIdxs                   = cp_idxs(gMRA.cp_nIdxs(idxs_kept)+1);                                                  % Update points kept
        gMRA.isaleaf                    = leafnodes(gMRA.cp)';
        gMRA.LeafNodes                  = find(gMRA.isaleaf);
        if gMRA.opts.verbose>0; fprintf('\n\t Construct cover tree on leaves...'); end
        gMRA.Timing.covertree_onleaves  = cputime;
        if isempty(gMRA.sorted_idxs)                                                                                            % Construct cover tree on leaves only
            opts.CoverTreeOpts.SamplesLeavesIdxs = gMRA.cp_nIdxs(gMRA.LeafNodes);
        else
            opts.CoverTreeOpts.SamplesLeavesIdxs = gMRA.sorted_idxs(gMRA.cp_nIdxs(gMRA.LeafNodes));
        end
        Tree_on_leaves                  = covertree_build( opts.CoverTreeOpts.Samples(:,opts.CoverTreeOpts.SamplesLeavesIdxs),opts.CoverTreeBuildOpts );           
        gMRA.Timing.covertree_onleaves  = cputime-gMRA.Timing.covertree_onleaves;
        if gMRA.opts.verbose>0; fprintf('done. (%.3f sec)',gMRA.Timing.covertree_onleaves); end
        if gMRA.opts.verbose>0; fprintf('\n\t Perform nearest neighbor searches on leaves...'); end
        gMRA.Timing.covertree_nnleaves  = cputime;                                                                          
        cmember                         = covertree_nnsearch( opts.CoverTreeOpts.Samples(:,opts.CoverTreeOpts.SamplesLeavesIdxs), Tree_on_leaves, X, 1 );     % Assigning points to leaves - MM:TBD:when Tree_on_leaves is empty!!
        gMRA.Timing.covertree_nnleaves  = cputime-gMRA.Timing.covertree_nnleaves;
        if gMRA.opts.verbose>0; fprintf('done. (%.3f sec)',gMRA.Timing.covertree_nnleaves); end
        gMRA.IniLabels                  = gMRA.LeafNodes(cmember)';
        gMRA.IniLabelsPtIdxs            = opts.CoverTreeOpts.SamplesLeavesIdxs(cmember)';
        gMRA.cp                         = gMRA.cp';
    otherwise
        error('Unknown PartitionType %s.',opts.PartitionType);
end

gMRA.Scales                 = Compute_scales(gMRA.cp);


if gMRA.opts.verbose,   fprintf('done. (%.3f sec)',toc); end

%% Initialize GWT structure
nAllNets        = length(gMRA.cp);
inputThreshold0 = gMRA.opts.threshold0;
J               = max(gMRA.Scales);
J0              = numel(inputThreshold0);
if J0 < J % too short, the use last element for finer scales
    inputThreshold0 = inputThreshold0*ones(1,J-J0+1);
end
% The following line converts threshold0 from per scale (1 by J vector) to per node (1 by nAllNets vector):
gMRA.opts.threshold0 = inputThreshold0(gMRA.Scales);                % node-dependent precision

% if gMRA.opts.precision<min(gMRA.opts.threshold0)      % separate thresholds for leaf nodes are provided
%     if gMRA.opts.pruning % minimal encoding cost pruning
%         gMRA.opts.threshold0(1:end) = gMRA.opts.precision;
%     else
%         gMRA.opts.threshold0(gMRA.LeafNodes) = gMRA.opts.precision;   % overwrite the leaf nodes with the specified precision
%     end
% end

% Initialize fields of gMRA
gMRA.Radii        = zeros(1,nAllNets);
gMRA.Sizes        = zeros(1,nAllNets);
gMRA.PointsInNet  = cell(1,nAllNets);
gMRA.Centers      = cell(1,nAllNets);
gMRA.ScalBasis    = cell(1,nAllNets);
gMRA.Sigmas       = cell(1,nAllNets);
gMRA.WavBases     = cell(1,nAllNets);
gMRA.WavConsts    = cell(1,nAllNets);
gMRA.WavSingVals  = cell(1,nAllNets);

if gMRA.opts.sparsifying,
    gMRA.Projections = cell(1,nAllNets); % coordinates of the projected data to each node
end

if gMRA.opts.splitting
    gMRA.WavDimsPsiCap = zeros(1,nAllNets); % intersection dimensions of the wavelet subspaces associated to the children of a specific node
end

% back up gMRA.cp
gMRA.cp_orig = gMRA.cp;

%% Construct geometric wavelets
% if gMRA.opts.parallel,
%     fprintf('\n');
%     delete(gcp('nocreate'));
%     gcp;
%     fprintf('\n');
% end;
if gMRA.opts.verbose,     fprintf('\n Constructing geometric wavelets...'); end
gMRA.Timing.GMRAconstruction = cputime;
if gMRA.opts.orthogonalizing
    gMRA = construct_orthogonalGeometricWavelets( X, gMRA );                   % orthogonal geometric wavelets
elseif gMRA.opts.pruning
    gMRA = construct_pruningGeometricWavelets( X, gMRA );                      % minimal encoding-cost pruning
else
    gMRA = construct_GMRA( X, gMRA );                                          % best approximations, no pruning
end
gMRA.Timing.GMRAconstruction = cputime - gMRA.Timing.GMRAconstruction;
if gMRA.opts.verbose,     fprintf('done. (%.3f sec) \n',gMRA.Timing.GMRAconstruction); end

%% Construct a cover tree on the centers of the leaves, which is going to be needed for computing FGWT's
if gMRA.opts.verbose,   fprintf('\n Constructing data structure for future FGWT''s...'); end
gMRA.Timing.Covertree_onleaves_centers      = cputime;
gMRA.opts.Covertree_on_leavesOpts           = gMRA.opts.CoverTreeBuildOpts;
%gMRA.opts.Covertree_on_leavesOpts.numlevels = gMRA.opts.CoverTreeTrimOpts.TrimValue;
[~, ~, ~, NNInfo] = nrsearch(cat(2,gMRA.Centers{gMRA.LeafNodes}), gMRA.Centers{gMRA.LeafNodes(1)}, 1, 0, struct('ReturnAsArrays',1,'CoverTreeBuildOpts',gMRA.opts.Covertree_on_leavesOpts));
gMRA.Timing.Covertree_onleaves_centers = cputime - gMRA.Timing.Covertree_onleaves_centers;
if isfield(NNInfo,'CoverTree'),     gMRA.TreeOnLeaves = NNInfo.CoverTree;   end
gMRA.Timing.PrecomputeFGWTdata = cputime;
gMRA = PrecomputeFGWTdata( gMRA );
gMRA.Timing.PrecomputeFGWTdata = cputime - gMRA.Timing.PrecomputeFGWTdata;
if gMRA.opts.verbose,     fprintf('done. (tree(%.3f sec), FGWT data(%.3f sec)) \n',gMRA.Timing.Covertree_onleaves_centers,gMRA.Timing.PrecomputeFGWTdata); end

%% Compute cost of dictionary
if gMRA.opts.verbose>0; fprintf('\n Computing cost of dictionary...'); end
gMRA.Timing.dictionarycosts = cputime;
gMRA.DictCosts              = ComputeDictionaryCost(gMRA);
gMRA.Timing.dictionarycosts = cputime - gMRA.Timing.dictionarycosts;
if gMRA.opts.verbose>0; fprintf('done.'); end

%%
gMRA.Timing.GW = cputime-gMRA.Timing.GW;

%% Shut down parallel pool
% delete(gcp('nocreate'));

gMRA.opts   = opts;

return;
