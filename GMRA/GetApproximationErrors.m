function [ApproxError,Quantiles] = GetApproximationErrors(X, gMRA, XGWT, Options)
%
% function [ApproxError,Quantiles] = GetApproximationErrors(X, gMRA, XGWT, Options)
%
%   Computes MSE across all scales
%
% IN:
%   X           : data set
%   gMRA        : GMRA data structure
%   XGWT        : transform of X using GMRA gMGRA
%   [cOpts]         : structure with the following fields:
%                       norm    : p-norm to be used to measure error at each scale, between cX and cProjections(:,:,j).
%                                 Can be any p\in(0,+\infty]. Note that norm=2 means the mean L^2 error, which is not the same
%                                 as the sqrt of MSE. It an be a vector. Default: 2.
%                       quantiles: list of quantiles of the errors. Default: [].
%
% OUT:
%   ApproxError.Absolute    : J by #Options.norm matrix whose entry (j,i) is the approximation error at scale j in norm Options.norm(i)
%   ApproxError.Relative    : J by #Options.norm matrix whose entry (j,i) is the approximation error at scale j in norm Options.norm(i)
%   Quantiles.Absolute      : J by #Options.norm by #Options.quantiles tensor whose entry (j,i,l) is the quantile Options.quantiles(l) of the error at scale j in norm Options.norm(i)
%   Quantiles.Relative      : J by #Options.norm by #Options.quantiles tensor whose entry (j,i,l) is the quantile Options.quantiles(l) of the error at scale j in norm Options.norm(i)
%

% (c) Mauro Maggioni, 2015

ApproxError = [];
Quantiles   = [];

% Parameters
if nargin<4,                      Options = [];                                 end
if ~isfield(Options,'norm'),      Options.norm = 2;                             end
if ~isfield(Options,'quantiles') || nargout<2, Options.quantiles = [];          end

% Memory allocation
ApproxError.Absolute = zeros(max(gMRA.Scales),length(Options.norm));
ApproxError.Relative = zeros(max(gMRA.Scales),length(Options.norm));
if ~isempty(Options.quantiles)
    Quantiles.Absolute      = zeros(max(gMRA.Scales),length(Options.norm),length(Options.quantiles));
    Quantiles.Relative      = zeros(max(gMRA.Scales),length(Options.norm),length(Options.quantiles));
end

% Pre-computations
lNormXk         = sqrt(sum(X.^2,1));
lNormXkisZero   = find(lNormXk==0);

% Go through the scales
for j = max(gMRA.Scales):-1:1
    Projectionj = IGWTScalej(gMRA,XGWT,j);                                                                                      % Compute data approximation at scale j
    
    lNormErrAbs = sqrt(sum((X-Projectionj).^2,1));
    
    for l = 1:length(Options.norm)
        lNormErrRel                 = lNormErrAbs./lNormXk;
        lNormErrRel(lNormXkisZero)  = lNormErrAbs(lNormXkisZero);

        switch Options.norm(l)                                                                                                  % Compute mean p-norm
            case inf
                ApproxError.Absolute(j,l)   = max(lNormErrAbs);                                                                 % Infinity norm
                ApproxError.Relative(j,l)   = max(lNormErrRel);
            otherwise
                lNormsp                     = lNormErrAbs.^(Options.norm(l));
                ApproxError.Absolute(j,l)   = (mean(lNormsp))^(1/Options.norm(l));                                              % Other p-norm
                lNormsp                     = lNormErrRel.^(Options.norm(l));
                ApproxError.Relative(j,l)   = (mean(lNormsp))^(1/Options.norm(l));                                              
        end
        if ~isempty(Options.quantiles)
            Quantiles.Absolute(j,l,:)        = quantile(lNormErrAbs,Options.quantiles);                                         % Quantiles
            Quantiles.Relative(j,l,:)        = quantile(lNormErrRel,Options.quantiles);
        end
    end
end

return