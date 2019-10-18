function [ApproxError,Quantiles] = ComputeAndDisplayApproxErr( X, gMRA, XGWT, ErrOpts, Sigma )

%
% function [ApproxError,Quantiles] = ComputeAndDisplayApproxErr( X, gMRA, XGWT, ErrOpts )
%
% Plots reconstruction error (w.r.t. to clean and noisy data) as a function of scale.
%
% USES:
%   GetApproximationErrors, DisplayApproximationErrors

if nargin<4 || isempty(ErrOpts),
    ErrOpts = struct('norm',[2,inf],'relative',true,'quantiles',[0.25,0.5,0.75]);                                       % Default option for approximation error calculation
end
if nargin<5 || isempty(Sigma)
    Sigma = 0;
end

% Compute approximation errors
[ApproxError,Quantiles] = GetApproximationErrors(X, gMRA, XGWT, ErrOpts );

return