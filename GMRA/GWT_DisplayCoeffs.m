function GWT_DisplayCoeffs( gMRA, XGWT )

%
% Displays tables of (1) magnitude of wavelet coefficients, (2) wavelet coefficients, (3) dimensions of the wavelet subspaces
%

J = max(gMRA.Scales);
N = sum(cellfun(@(B)size(B,2), XGWT.CelScalCoeffs(:,1) ));

%% Form the matrix of wavelet coefficients
if ~isfield(XGWT,'MatWavCoeffs'), 
    [XGWT.MatWavCoeffs, XGWT.maxWavDims, XGWT.MatWavDims] = XGWT2mat( gMRA, XGWT );                                             % Reshape FGWT into matrix if desired
end
WavCoeffs   = XGWT.MatWavCoeffs;
WavCoeffs2  = WavCoeffs.^2; 

WavCoeffMags = zeros(N,J);
for j = 1:J
   WavCoeffMags(:,j) = sqrt(sum(WavCoeffs2(sum(XGWT.maxWavDims(1:j-1))+1:sum(XGWT.maxWavDims(1:j)),:),1));
end

%% Display wavelet coefficient magniture and dimension of wavelet subspaces
DrawGWTWavSubFcn( gMRA, log10(WavCoeffMags)/log10(1/gMRA.TreeOnLeaves.theta), 'Magnitude of Wavelet Coefficients (log_{1/\theta} scale)' );
DrawGWTWavSubFcn( gMRA, XGWT.MatWavDims, 'Dimension of wavelet subspaces' );

%% Display Coefficients against scales
MeanWavCoeffMags = mean(WavCoeffMags, 1)';

% delta = zeros(1,J);
% for j = 1:J
%     delta(j) = mean(gMRA.Radii(gMRA.Scales == j));
% end

delta = gMRA.Radii_AvePerScale(unique(gMRA.Scales));

figure; 
goodidxs = find(MeanWavCoeffMags>2*eps*max(MeanWavCoeffMags));
if J>3 
    PlotPolyFit( -log10(delta(2:max(goodidxs))')/log10(1/gMRA.TreeOnLeaves.theta), log10(MeanWavCoeffMags(2:max(goodidxs)))/log10(1/gMRA.TreeOnLeaves.theta) ); 
elseif J>1
    PlotPolyFit( -log10(delta')/log10(1/gMRA.TreeOnLeaves.theta), log10(MeanWavCoeffMags)/log10(1/gMRA.TreeOnLeaves.theta) );
end
plot(-log10(delta)/log10(1/gMRA.TreeOnLeaves.theta), log10(MeanWavCoeffMags)/log10(1/gMRA.TreeOnLeaves.theta), 'o-','LineWidth',1,'Color',[0,0,0]); 
title('Wavelet coefficients vs. scale (log_{1/\theta} scale)', 'fontSize', 12)
xlabel('scale', 'fontSize', 12); ylabel('coefficient', 'fontSize', 12)
grid on

axis tight;

return


function DrawGWTCoeffs( gMRA, Coeffs, maxWavDims, Title )

if nargin<4, Title =''; end

J = max(gMRA.Scales);

%% Show the matrix of all wavelet coefficients
figure; imagesc(Coeffs); 
cmap=map2;
colormap(cmap); 
title(Title, 'fontSize', 12);
xlabel('points', 'fontSize', 12); ylabel('scales', 'fontSize', 12);
linewdth = 1;
col = 0.8*ones(1,3);

for j = 1:J
    nets = sort([find(gMRA.Scales == j);gMRA.LeafNodes(gMRA.Scales(gMRA.LeafNodes)<j)], 'ascend');
    lens  = cumsum(gMRA.Sizes(nets)); 
    ln = line([0 lens(end)], [sum(maxWavDims(1:j-1))+0.5 sum(maxWavDims(1:j-1))+0.5]); set(ln,'color',col,'LineWidth',linewdth);
    for i = 1:length(lens)-1
        ln = line([lens(i) lens(i)], [sum(maxWavDims(1:j-1))+0.5 sum(maxWavDims(1:J))+0.5]); set(ln,'color',col,'LineWidth',linewdth);
    end
end
balanceColor
colorbar

return;

