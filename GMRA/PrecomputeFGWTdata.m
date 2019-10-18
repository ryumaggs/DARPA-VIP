function gMRA = PrecomputeFGWTdata( gMRA )

% 
% function gMRA = PrecomputeFGWTdata( gMRA )
% Pre-computes some change-of-basis matrices to be used for FGWT

% List nodes and their parents in top-bottom fashion
[~,nodeList] = sort(gMRA.Scales);
parentList   = gMRA.cp(nodeList);

%
gMRA.ScalCenter                 = cell(length(nodeList),1);
gMRA.ScalBasisChange            = cell(length(nodeList),1);
gMRA.ScalWavBasisChange         = cell(length(nodeList),1);
gMRA.ScalCenterConsts           = cell(length(nodeList),1);
gMRA.ScalBasisCoarserCenter     = cell(length(nodeList),1);

gMRA.WavConstsCoeffs            = cell(length(nodeList),1); 

for k = 1:length(nodeList)
    if ~isempty(gMRA.ScalBasis{nodeList(k)})
        gMRA.ScalCenter{nodeList(k)}                = gMRA.ScalBasis{nodeList(k)}*gMRA.Centers{nodeList(k)};
    end
    
    if parentList(k)>0 && ~isempty(gMRA.ScalBasis{parentList(k)})
        if ~isempty(gMRA.ScalBasis{nodeList(k)}),
            gMRA.ScalBasisChange{nodeList(k)}       = gMRA.ScalBasis{parentList(k)}*gMRA.ScalBasis{nodeList(k)}';
        end
        if ~isempty(gMRA.WavBases{nodeList(k)})
            gMRA.ScalWavBasisChange{nodeList(k)}    = gMRA.ScalBasis{parentList(k)}*gMRA.WavBases{nodeList(k)}';
        end
        if ~isempty(gMRA.WavConsts{nodeList(k)})
            gMRA.ScalCenterConsts{nodeList(k)}      = gMRA.ScalBasis{parentList(k)}*(gMRA.Centers{nodeList(k)}+gMRA.WavConsts{nodeList(k)});
        else
            gMRA.ScalCenterConsts{nodeList(k)}      = gMRA.ScalBasis{parentList(k)}*gMRA.Centers{nodeList(k)};
        end
        gMRA.ScalBasisCoarserCenter{nodeList(k)}    = gMRA.ScalBasis{parentList(k)}*gMRA.Centers{nodeList(k)};
        
        if ~isempty(gMRA.WavBases{nodeList(k)})
            gMRA.WavConstsCoeffs{nodeList(k)}       = gMRA.WavBases{nodeList(k)}*gMRA.WavConsts{nodeList(k)};
        end
    end
end

return
