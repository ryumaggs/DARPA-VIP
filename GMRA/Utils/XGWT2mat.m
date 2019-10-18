function [MatWavCoeffs, maxWavDims, MatWavDims] = XGWT2mat( gMRA, XGWT )

%
% function [MatWavCoeffs, maxWavDims, MatWavDims] = XGWT2mat( gMRA, XGWT )
%
% IN:
%   XGWT    : FGWT of a data set with N points
%
% OUT:
%   MatWavCoeffs    : A sum(maxWavDims)xN matrix of wavelet coefficients at all scales (from coarse to fine)
%   maxWavDims      : J vector of maximum dimension of wavelet subspaces at each scale
%   MatWavDims      : A JxN matrix of dimensions of wavelet subspaces
%

% (c) Mauro Maggioni, 2015, modified from initial version by G. Chen.

J           = size(XGWT.CelWavCoeffs,2);                                                                                        %% Determine various parameters
szBlocks    = cellfun(@(B)size(B,2), XGWT.CelWavCoeffs(:,1) );                                                                  % Compute sizes of all the leaves
cumSzBlocks = [1;1+cumsum(szBlocks)];
N           = cumSzBlocks(end)-1;

MatWavCoeffs= cell(J,1);                                                                                                        % Memory allocation
maxWavDims  = zeros(1,J);                                                                                                        
MatWavDims  = zeros(N,J);                                                                                                        

for j = 1:J                                                                                                                     %% Go through all scales    
    BlockWavDims_j  = cellfun(@(B)size(B,1), XGWT.CelWavCoeffs(:,j));                                                           % Each block is a leaf node
    maxWavDims(j)   = max(BlockWavDims_j);                                                                                      % Maximum dimension of wavelet subspaces at this scale
    MatWavCoeffs{j} = zeros(maxWavDims(j),N);    
    positiveWavDims = (find(BlockWavDims_j>0))';
    
    for n = positiveWavDims                                                                                                     % Go through non-empty wavelet subspaces
        MatWavCoeffs{j}(1:BlockWavDims_j(n),cumSzBlocks(n):cumSzBlocks(n+1)-1)  = XGWT.CelWavCoeffs{n,j};                       % Copy wavelet coefficients into set of coefficients at scale j
        MatWavDims(cumSzBlocks(n):cumSzBlocks(n+1)-1,j)                         = BlockWavDims_j(n);                            % Save wavelet subspace dimension inot matrix of wavelet dimensions
    end    
end

MatWavCoeffs = cat(1,MatWavCoeffs{:});                                                                                          % Group all wavelet coefficients at all scales into matrix

return
