function [Net,Weights,Density] = ConstructNet( X, r, opts )

%
% function Net=ConstructNet( X, r )
%
% IN:
%   X       : D by N matrix of data points
%   r       : radius for net
%   [opts]  : structure of options:
%               [weightfcn] : weight function to apply to distances of points. For example @exp(-dist^2/sigma)
%
% OUT:
%   Net     : indices (into columns of X) of the net points
%
% COMMENTS:
%   MM: should scale roughly linearly in N, and linearly in the size of the net. I have only done small tests (up to 10^6 points).
%
% EXAMPLE
%   X = randn(2,100000); 
%   [Net,Weights,Density]=ConstructNet(X,0.2); 1/size(X,2)*X*X',(1/sum(Weights)*X(:,Net)*diag(single(Weights))*(X(:,Net))'),
%   [Net,Weights,Density]=ConstructNet(X,0.2,struct('weightfcn',@(x) exp(-x.^2/0.1)));
%   figure;plot(X(1,:),X(2,:),'.');hold on;scatter(X(1,Net),X(2,Net),20,Weights,'filled');colorbar; figure;plot(X(1,:),X(2,:),'.');hold on;scatter(X(1,Net),X(2,Net),20,Density,'filled');colorbar;
%   
%
% (c) Mauro Maggioni
%
FRACTION_FOR_RECURSION = 9/10;

if nargin<3,                    opts = struct();                                        end
if nargin<2 || isempty(r),      r    = EstimateRforNet(X, opts);                        end
if isfield(opts,'weightfcn'),   weightFcn = opts.weightfcn;     else weightFcn = [];    end

N           = size(X,2);
Net         = zeros(N,1,'uint32');
Weights     = zeros(N,1,'uint32');
activeIdxs  = ones(1,N)>0;
nactiveIdxs = N;

X           = X';
curidx      = 1;

rndperm     = randperm(N);
for k = 1:N
    found = false;
    for i = k:N                                                                                                                 % Find a new net point
        Net(curidx) = rndperm(i);
        if activeIdxs(Net(curidx)), found = true; break; end;
    end             
    if ~found,  continue;   end
    dists                               = pdist2(X(Net(curidx),:),X);                                                           % Compute distance between net point and all other points
    dists_withinradiusIdxs              = dists<r;                                                                              % Gets points within specified radius
    dists_withinradiusIdxs              = dists_withinradiusIdxs & activeIdxs;
    if isempty(weightFcn)
        Weights(curidx)                 = sum(dists_withinradiusIdxs);                                                          % Weight is simply the number of active points within specified radius
    else
        Weights(curidx)                 = sum(weightFcn(dists(dists_withinradiusIdxs)));                                        % Weight is the sum of weights computed by weightfcn applied to distances to active points within specified radius
    end
    activeIdxs(dists_withinradiusIdxs)  = false;
    nactiveIdxs                         = nactiveIdxs - Weights(curidx);
    curidx                              = curidx+1; 
    if nactiveIdxs < N*FRACTION_FOR_RECURSION                                                                                   % Recurse on smaller subset
        [newNet,newWeights]                 = ConstructNet( X(activeIdxs,:)', r );
        activeIdxs_idxs                     = find(activeIdxs);
        Net(curidx:curidx+length(newNet)-1) = activeIdxs_idxs(newNet);
        Weights(curidx:curidx+length(newNet)-1) = newWeights;
        curidx                              = curidx + length(newNet);
        break
    end
end

Net(curidx:end) = [];
Weights(curidx:end) = [];

if nargout>2                                                                                                                    % Compute weights if requested
    dists                               = pdist2(X(Net,:),X);                                                                   % Compute distance between net points and all other points
    Density                             = sum(dists<r/2,2);
end

return


function r = EstimateRforNet( X, opts )

if nargin<2,                            opts = struct();                                end
if ~isfield(opts,'downsamplingFactor'), opts.downsamplingFactor = 10;                   end

N   = size(X,2);
N0  = ceil(100*log(N));

X   = X';

rndperm     = randperm(N);

dists   = pdist2(X(rndperm(1:N0),:),X);
for k = 1:size(dists,1)
    dists(k,:)   = sort(dists(k,:),'ascend');
end

r       = median(dists(:,min([opts.downsamplingFactor,size(dists,2)])));

return