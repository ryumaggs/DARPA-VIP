function [ApproxError, Quantiles, X_Approx] = GetApproximationError(gMRA, X, XGWT, Partition, Options )

%
% function GetApproximationError(gMRA, X, XGWT, Partition, Options )
%
%   Computes approximation error corresponding to a partition
%
% IN:
%   gMRA        : GMRA data structure
%   X           : data set
%   XGWT        : transform of X using GMRA gMGRA
%   Partition   : partition, i.e. vector of indices of gMRA nodes
%   [Options]   : structure with the following fields:
%                       norm    : p-norm to be used to measure error at each scale, between cX and cProjections(:,:,j).
%                                 Can be any p\in(0,+\infty], and equal to (\frac 1N \sum_{i=1}^N ||x_i-\hat x_i||_2^p)^{\frac 1p}
%                       quantiles: list of quantiles of the errors. Default: [].
%
% OUT:
%   ApproxError.Absolute    : J by #Options.norm matrix whose entry (j,i) is the approximation error at scale j in norm Options.norm(i)
%   ApproxError.Relative    : J by #Options.norm matrix whose entry (j,i) is the approximation error at scale j in norm Options.norm(i)
%   Quantiles.Absolute      : J by #Options.norm by #Options.quantiles tensor whose entry (j,i,l) is the quantile Options.quantiles(l) of the error at scale j in norm Options.norm(i)
%   Quantiles.Relative      : J by #Options.norm by #Options.quantiles tensor whose entry (j,i,l) is the quantile Options.quantiles(l) of the error at scale j in norm Options.norm(i)
%   X_Approx                : approximation of X on Partition
%

% (c) Mauro Maggioni, 2015

ApproxError = [];
Quantiles   = [];

% Parameters
if nargin<4,                      Options = [];                                 end
if ~isfield(Options,'norm'),      Options.norm = 2;                             end
if ~isfield(Options,'quantiles') || nargout<2, Options.quantiles = [];          end

% Memory allocation
ApproxError.Absolute = zeros(1,length(Options.norm));
ApproxError.Relative = zeros(1,length(Options.norm));
if ~isempty(Options.quantiles)
    Quantiles.Absolute      = zeros(length(Options.norm),length(Options.quantiles));
    Quantiles.Relative      = zeros(length(Options.norm),length(Options.quantiles));
end

X_Approx        = GetApproximationOnPartition(gMRA,Partition,XGWT,struct('type',Options.type));                                                             % Get approximation on partition

lNormXk         = sqrt(sum(X.^2,1));                                                                                            % Pre-computations
lNormXkisZero   = find(lNormXk==0);

lNormErrAbs = sqrt(sum((X-X_Approx).^2,1));                                                                                     % Approximation error...

for l = 1:length(Options.norm)                                                                                                  %...in different norms
    lNormErrRel                         = lNormErrAbs./lNormXk;
    lNormErrRel(lNormXkisZero)          = lNormErrAbs(lNormXkisZero);
    
    switch Options.norm(l)                                                                                                      % Compute mean p-norm
        case inf
            ApproxError.Absolute(l)     = max(lNormErrAbs);                                                                     % Infinity norm
            ApproxError.Relative(l)     = max(lNormErrRel);
        otherwise
            lNormsp                     = lNormErrAbs.^(Options.norm(l));
            ApproxError.Absolute(l)     = (mean(lNormsp))^(1/Options.norm(l));                                                  % Other p-norm
            lNormsp                     = lNormErrRel.^(Options.norm(l));
            ApproxError.Relative(l)     = (mean(lNormsp))^(1/Options.norm(l));
    end
    if ~isempty(Options.quantiles)
        Quantiles.Absolute(l,:)         = quantile(lNormErrAbs,Options.quantiles);                                              % Quantiles
        Quantiles.Relative(l,:)         = quantile(lNormErrRel,Options.quantiles);
    end
end

return