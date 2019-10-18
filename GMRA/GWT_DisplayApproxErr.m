function [ApproxError,Quantiles] = ComputeAndDisplayApproxErr( X, gMRA, XGWT, ErrOpts )

%
% function [ApproxError,Quantiles] = ComputeAndDisplayApproxErr( X, gMRA, XGWT, ErrOpts )
%
% Plots reconstruction error (w.r.t. to clean and noisy data) as a function of scale.
%
% USES:
%   GetApproximationErrors, DisplayApproximationErrors

if nargin<4,
    ErrOpts = struct('norm',[2,inf],'relative',true,'quantiles',[0.1,0.25,0.5,0.75,0.9]);                                       % Default option for approximation error calculation
end

% Compute approximation errors
[ApproxError,Quantiles] = GetApproximationErrors(X, gMRA, XGWT, ErrOpts );

% Display approximation error
DisplayApproximationErrors( gMRA,ApproxError,Quantiles);

return