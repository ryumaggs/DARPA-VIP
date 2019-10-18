function X_Approx = FAGWT_Approximation(gMRA,Partition,XGWT)

%% Perform adaptive GMRA
% Step 2:   Build piecewise linear approximation of the test samples on the adaptive partition
% Partition: partition
% XGWT_testtest: Scaling coefficients for the test samples

n  = size(XGWT.DeltaPoint,1);
if ~isfield(gMRA.opts,'Proj') || isempty(gMRA.opts.Proj)
    X_Approx    = zeros(max(cellfun(@(x) size(x,2),gMRA.ScalBasis)),sum(XGWT.leafNodeSizes));                                   %sum(cellfun(@(x) size(x,2),XGWT.CelScalCoeffs(:,1))));
    ProjT       = [];
else
    X_Approx    = zeros(size(gMRA.opts.Proj,2),sum(XGWT.leafNodeSizes));                                                        %sum(cellfun(@(x) size(x,2),XGWT.CelScalCoeffs(:,1))));
    ProjT       = gMRA.opts.Proj';
end


% Build approximation
nLeafNodes              = numel(gMRA.LeafNodes);

for i = 1: nLeafNodes
    iFineNet            = Partition(i);
    j_loc               = gMRA.Scales(Partition(i));
    if ~isempty(XGWT.CelScalCoeffs{i,j_loc})
        proj_i = bsxfun(@plus, (gMRA.ScalBasis{iFineNet})'*XGWT.CelScalCoeffs{i,j_loc}, gMRA.Centers{iFineNet});
        if ~isempty(ProjT)
            proj_i = ProjT*proj_i;
        end        
        if isfield(gMRA.opts,'Mean') && ~isempty(gMRA.opts.Mean)
            proj_i = bsxfun(@plus,proj_i,gMRA.opts.Mean);
        end        
        X_Approx(:,XGWT.PointsInNet{gMRA.LeafNodes(i)}) = proj_i;
    end
end

return
