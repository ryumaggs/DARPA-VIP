function X_Approx = GetApproximationOnPartition(gMRA,Partition,XGWT,opts)

%
% function X_Approx = GetApproximationOnPartition(gMRA,Partition,XGWT,opts)
%
% IN:
%   gMRA      : GMRA sturcture
%   Partition : partition
%   XGWT      : Scaling coefficients for the samples
% OUT:
%   X_Approx  : DxN matrix representing the approximations to the points
%   opts      : opts.type == GMRA    GMRA approximation, default
%               opts.type == Center  approximation by the center
%     
%
% (c) Mauro Maggioni, Duke University
%

if nargin == 3
    opts.type = 'GMRA';
end

n  = size(XGWT.DeltaPoint,1);
if ~isfield(gMRA.opts,'Proj') || isempty(gMRA.opts.Proj)
    X_Approx    = zeros(gMRA.opts.EmbedDim,sum(XGWT.leafNodeSizes));                                   %sum(cellfun(@(x) size(x,2),XGWT.CelScalCoeffs(:,1))));
    ProjT       = [];
else
    X_Approx    = zeros(size(gMRA.opts.Proj,2),sum(XGWT.leafNodeSizes));                                                        %sum(cellfun(@(x) size(x,2),XGWT.CelScalCoeffs(:,1))));
    ProjT       = gMRA.opts.Proj';
end


% Build approximation

for i = 1:length(Partition)
    iFineNet            = Partition(i);
    j_loc               = gMRA.Scales(iFineNet);
    if length(j_loc)>2, warning('ops!'), keyboard; end
    if strcmp(opts.type,'GMRA') && ~isempty(XGWT.CelScalCoeffs{i,j_loc})
        proj_i = bsxfun(@plus, gMRA.ScalBasis{iFineNet}'*XGWT.CelScalCoeffs{i,j_loc}, gMRA.Centers{iFineNet});
    elseif strcmp(opts.type,'Center') || isempty(XGWT.CelScalCoeffs{i,j_loc})
        proj_i = repmat(gMRA.Centers{iFineNet},1,length(XGWT.PointsInNet{gMRA.LeafNodes(i)}));
    end
    if ~isempty(ProjT)
        proj_i = ProjT*proj_i;
    end
    if isfield(gMRA.opts,'Mean') && ~isempty(gMRA.opts.Mean)
        proj_i = bsxfun(@plus,proj_i,gMRA.opts.Mean);
    end
    X_Approx(:,XGWT.PointsInNet{gMRA.LeafNodes(i)}) = proj_i;
end

return
