function nodeInfo = local_SVD_analysis( X, gMRA,node,ptsInNode_idxs )
%
% function nodeInfo = local_SVD_analysis( X, gMRA,node,ptsInNode_idxs )
%
% IN:
% X         : D by N data matrix
% gMRA      : GMRA data structure
% node      : node to analyse
% ptsInNode_idxs : indices into columns of X for which points are assigned to this node
%
% OUT:
% nodeInfo  :   structure with the following fields:
%               Size    : number of points in node
%               Center  : mean of the points in node
%               Radius  : radius of points in node
%               Sigmas  : singular values (->nodeInfo.Sigmas
%               ScalBasis: singular vectors (->nodeInfo.ScalBasis
%               epsEncodingCost : cost of encoding points in node (->nodeInfo.epsEncodingCost
%               Projections : projection of points on scaling functions (->gMRA.Projections{node} (only if sparsifying

MAX_PTS_FOR_MEAN    = Inf;

% This function performs local svd analysis at the input node.
% In particular, it computes local mean and basis vectors, as well as the
% singular values.

nodeInfo.Size   = numel(ptsInNode_idxs);
Xlocal          = X(:,ptsInNode_idxs);
if isempty(Xlocal),
    %fprintf('\n WARNING: local_SVD_analysis: empty node %d',node);
    nodeInfo.Center     = [];
    nodeInfo.Radius     = 0;
    nodeInfo.Sigmas     = [];
    nodeInfo.ScalBasis  = [];
    if gMRA.opts.pruning, nodeInfo.epsEncodingCost = 0; end
    return;
end

if MAX_PTS_FOR_MEAN==Inf
    nodeInfo.Center     = mean(Xlocal,2);
else
    lPtsInNode_rand_idxs= randperm(length(ptsInNode_idxs));
    lPtsInNode_rand_idxs= lPtsInNode_rand_idxs(1:min([MAX_PTS_FOR_MEAN,length(lPtsInNode_rand_idxs)]));
    lPtsInNode_rand     = Xlocal(:,lPtsInNode_rand_idxs);
    nodeInfo.Center     = mean(lPtsInNode_rand,2);
end

Y = bsxfun(@minus,Xlocal,nodeInfo.Center);                                                                                      % Centered data in the current node

nodeInfo.Radius    = sqrt(max(sum(Y.^2,1)));

%% Compute local SVD
if gMRA.opts.ManifoldDimension==0 || gMRA.isaleaf(node)                                                                                              % Local dimension is not fixed, but based on local singular value decay                                                                  %% Local dimension is not fixed, but based on local singular value decay    
    [V,S,~]             = randPCA(Y,min([min(size(Y)),gMRA.opts.MaxDim]));                                                      % Use randomized PCA
    remEnergy           = sum(sum(Y.^2))-sum(diag(S).^2);
    nodeInfo.Sigmas     = ([diag(S); sqrt(remEnergy)]) /sqrt(nodeInfo.Size);
    if gMRA.opts.pruning %gMRA.isaleaf(node) || gMRA.opts.pruning
        errorType = gMRA.opts.errorType;
    else
        errorType = 'relative';
    end
    if ~gMRA.isaleaf(node)
        reqDim = min(numel(diag(S)), mindim(nodeInfo.Sigmas, errorType, gMRA.opts.threshold0(node)));
    else
        reqDim = min(numel(diag(S)), mindim(nodeInfo.Sigmas, errorType, gMRA.opts.precision));
    end
    nodeInfo.ScalBasis = V(:,1:reqDim);
else
    % Manifold dimension is given
    [V,S,~]             = randPCA(Y,min([gMRA.opts.ManifoldDimension,min(size(Y))]));                                           % Use randomized PCA
    nodeInfo.Sigmas     = diag(S)/sqrt(nodeInfo.Size);
    if size(V,2)<gMRA.opts.ManifoldDimension,
        V = [V,zeros(gMRA.opts.AmbientDimension,gMRA.opts.ManifoldDimension-size(V,2))];
    end
    nodeInfo.ScalBasis = V(:,1:min(gMRA.opts.ManifoldDimension,length(find(nodeInfo.Sigmas))));
end

%% Pruning
if gMRA.opts.pruning, % minimal encoding cost pruning
    nodeInfo.epsEncodingCost = (nodeInfo.Size+gMRA.opts.AmbientDimension) * reqDim + gMRA.opts.AmbientDimension;
end

%% Sparsify the local dictionary if requested
if gMRA.opts.sparsifying,
    %     gMRA.Projections{node} = Y*nodeInfo.ScalBasis*(nodeInfo.ScalBasis)';
    if gMRA.isaleaf(node) || (~gMRA.isaleaf(node) && gMRA.opts.addTangentialCorrections)
        gMRA.Projections{node} = nodeInfo.ScalBasis*nodeInfo.ScalBasis'*Y + repmat(nodeInfo.Center, 1,nodeInfo.Size);
    else %~isempty(children) && ~gMRA.opts.addTangentialCorrections
        children = (gMRA.cp==node);
        childrenProj = cat(2, gMRA.Projections{children});
        if ~isempty(childrenProj)
            gMRA.Projections{node} = nodeInfo.ScalBasis*(nodeInfo.ScalBasis)'*(childrenProj-repmat(nodeInfo.Center, 1,nodeInfo.Size)) + repmat(nodeInfo.Center, 1,nodeInfo.Size);
        end
        %computeWaveletCoeffcients(cat(2, gMRA.Projections{children})-repmat(nodeInfo.Center, 1,nodeInfo.Size), nodeInfo.ScalBasis, gMRA.opts.sparsifying)
    end
    %     if (nodeInfo.Size>gMRA.opts.ManifoldDimension) && (size(Y,1)>10),
    %         nodeInfo.ScalBasis = ksvd(struct('data', Y', 'Tdata', gMRA.opts.ManifoldDimension, 'initdict', nodeInfo.ScalBasis, 'iternum', 10));   %KSVD
    %         param.K=min([ceil(size(Y,1)/2),2*size(nodeInfo.ScalBasis,2)]);  % learns a dictionary with 100 elements
    %         param.lambda= 5;
    %         param.numThreads=-1; % number of threads
    %         param.approx=4.0;   %MM: no idea what this is for
    %         param.iter = -5; % let's wait only 20 seconds
    %         param.mode = 0;
    %         param.D = nodeInfo.ScalBasis;
    %         D = mexTrainDL(Y',param);
    %         param.mode = 1;
    %         param.lambda = norm(Y'-nodeInfo.ScalBasis*nodeInfo.ScalBasis'*Y');
    %         newcoeffs = mexLasso(Y',D,param);
    %         newcoeffs2 = mexLasso(Y',nodeInfo.ScalBasis,param);
    %         nodeInfo.ScalBasis = D;
    %    end
end

return
