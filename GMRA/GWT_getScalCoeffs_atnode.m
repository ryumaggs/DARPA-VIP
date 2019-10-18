function [ScalCoeffs,tangentialCorrections,k_idxs,j_idxs] = GWT_getScalCoeffs_atnode( gMRA, XGWT, cp_idx )

% Find the scaling coefficients of points belonging to each node at scale j
j_idxs                  = gMRA.Scales(cp_idx);
k_idxs                  = find(XGWT.Cel_cpidx(:,j_idxs)==cp_idx);
ScalCoeffs              = [];
tangentialCorrections   = [];

for i = 1:length(k_idxs),
    ScalCoeffs            = [ScalCoeffs;XGWT.CelScalCoeffs{k_idxs(i),j_idxs}'];
    if isfield(XGWT,'CelTangCoeffs') && nargout>1
        tangentialCorrections = [tangentialCorrections,XGWT.CelTangCoeffs{k_idxs(i),j_idxs}'];
    end
end

return