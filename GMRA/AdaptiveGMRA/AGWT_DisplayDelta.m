function AGWT_DisplayDelta(XGWT,UniformScales)
%
% function to display the approximation variations \|P_j x - P_{j+1} x\| at each point
% 
% Usage: AGWT_DisplayDelta(XGWT,UniformScales)
%
% Input: XGWT: FGWT of data
%        UniformScales: vector of length equal to the maximal scale, with
%        the jth entry being the scale at scale j
%
% (c) Wenjing Liao, Mauro Maggioni

if UniformScales(end) == inf
    UniformScales(end) = UniformScales(end-1);
end
figure
imagesc(1:size(XGWT.DeltaPoint,1),1:length(UniformScales),log10((XGWT.DeltaPoint')))
colorbar
title('log_{10}||P_{j+1}x-P_j x||^2 ')
xlabel('Point index')
ylabel('Scale')
caxis([-10,2])


% %% Display Deltajk at each cell
% figure
% imagesc(1:size(XGWT.DeltaCell,1),1:length(UniformScales),log10(XGWT.DeltaCell'))
% colorbar
% title('Cell-wise log_{10}||P_{j+1}-P_j||^2 ')
% xlabel('Cell index')
% ylabel('Scale')
% caxis([-10,2])


