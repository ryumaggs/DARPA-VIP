function [Projectionj,Partitionj] = IGWTScalej(gMRA,XGWT,j)

%
% function Projectionj = IGWTScalej(gMRA,XGWT,j)
%
% Compute the projection at Scale j
% Input: gMRA: GMRA structure
%        XGWT: Geometric wavelet transform data
%           j: Scale
% Output: Projection at scale j, partition at scale j (a vector of node indices into the gMRA nodes)
%
 
% (c) Mauro Maggioni, 2015
 
% Set up & memory allocation
nLeafNodes     = numel(gMRA.LeafNodes);

if ~isfield(gMRA.opts,'Proj') || isempty(gMRA.opts.Proj)
    Projectionj    = zeros(max(cellfun(@(x) size(x,2),gMRA.ScalBasis)),sum(cellfun(@(x) size(x,2),XGWT.CelScalCoeffs(:,1))));
    ProjT          = [];
else
    Projectionj    = zeros(size(gMRA.opts.Proj,2),sum(cellfun(@(x) size(x,2),XGWT.CelScalCoeffs(:,1))));
    ProjT          = gMRA.opts.Proj';
end
 
Partitionj  = zeros(nLeafNodes,1);

% Go through the leaf nodes
for i = 1:nLeafNodes
    if sum(XGWT.leafNodeIdxs==i)>0
        if isempty(XGWT.CelScalCoeffs{i,j})
            idx = find(XGWT.Cel_cpidx(i,1:j),1,'last');
        else
            idx = j;
        end
        iFineNet                                            = XGWT.Cel_cpidx(i,idx);
        proj_j                                              = bsxfun(@plus, (gMRA.ScalBasis{iFineNet})'*XGWT.CelScalCoeffs{i,idx}, gMRA.Centers{iFineNet});

        % Add mean and re-embed in high dimensions if needed
        if ~isempty(ProjT)
            proj_j = ProjT*proj_j;
        end
        if isfield(gMRA.opts,'Mean') && ~isempty(gMRA.opts.Mean)
            proj_j = bsxfun(@plus,proj_j,gMRA.opts.Mean);
        end
        
        Projectionj(:,XGWT.PointsInNet{gMRA.LeafNodes(i)})  = proj_j;
        Partitionj(i)                                       = iFineNet;
    end
end

Partitionj = unique(Partitionj);
 
return