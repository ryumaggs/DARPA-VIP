function Scores = ComputeRegressionScores( gMRAReg, X_test, XGWT, ScalCoeffsExt )

Npartition = length(XGWT.PointsInNet);

Scores = cell(length(gMRAReg.Reg),Npartition);

for deg = 1:length(gMRAReg.Reg)
    for p = 1:Npartition        
        if ~gMRAReg.RegressionOpts.noProj
            X_local     = ScalCoeffsExt{p}; %gMRAReg.gMRA.ScalBasis{p}*bsxfun(@minus,X_local,gMRAReg.gMRA.Centers{p});
            if isempty(X_local)
                X_local = zeros(1,length(XGWT.PointsInNet{p}));
            end
        else
            X_local = X_test(:,XGWT.PointsInNet{p});
        end
        poly = gMRAReg.Reg{deg}{p};
        
        if isfield(poly,'Lim') && poly.Lim>1
            for k = size( X_local,2 ):-1:1
                tmpScores               = repmat(X_local(:,k)',[length(poly.PowerMatrix) 1]).^poly.PowerMatrix;
                Scores{deg,p}(:,k)      = prod(tmpScores,2);
            end
        end
    end
end

return
