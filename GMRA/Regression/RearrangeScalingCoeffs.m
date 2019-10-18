function ScalCoeffsExt = RearrangeScalingCoeffs( gMRA, XGWT )

ScalCoeffsExt = cell(length(gMRA.cp),1);

for k = 1:length(gMRA.cp)
    ScalCoeffsExt{k} = GWT_getScalCoeffs_atnode( gMRA, XGWT, k )';                                    
end

return