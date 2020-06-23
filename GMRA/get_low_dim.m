function [low_dim] = get_low_dim(gMRA, data)
    low_dim = {};
    for i = 1 : size(data,2)
        xgwt = FGWT(gMRA, data(:, i));
        leafIndx = xgwt.leafNodeIdxs;
        low_dim{i} = xgwt.CelWavCoeffs(leafIndx,:);
        if rem(i, 100) == 0
            fprintf("checkpoint: %.1f\n", i)
        end
    end