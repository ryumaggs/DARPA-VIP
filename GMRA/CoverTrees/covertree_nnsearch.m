function [idxs,dists,opts] = covertree_nnsearch( X, CoverTree, Xq, kNN, opts )

%
% function [idxs,dists] = covertree_nnsearch( X, CoverTree, Xq, kNN, opts )
%
% Wrapper for cover tree findwithinMEX function
%
% IN:
%   X               : D by N matrix of N points in R^D
%   CoverTree       : covertree on cX, constructed by covertree
%   Xq              : D by M matrix M query points in R^D
%   kNN             : how many nearest neighbors to find
%   [opts]          : structure of options.
%                       [NTHREADS]      : number of threads for search. Default: number of cores.
%
% OUT:
%   idxs            : kNN by M array with the k-th column being the indices (into columns of cX) of the kNN nearest neighbors
%   dists           : kNN by M array of distances of the kNN nearest neighbors.
%   The above outputs are sorted in increasing distance.
%
%
% (c) Copyright Duke University, 2013
% Mauro Maggioni
% mauro@math.duke.edu
%
persistent NTHREADS

MIN_N_FOR_MULTITHREADED  = 1024;

if nargin<5,                    opts = struct();                                    end
if isunix && ~ismac
    NTHREADS = int32(0);
else    
    if size(X,2)>MIN_N_FOR_MULTITHREADED
        if isfield(opts,'NTHREADS'),        NTHREADS = int32(opts.NTHREADS);
        else                                NTHREADS = int32(feature('numcores'));      end
    else
        if isfield(opts,'NTHREADS'),        NTHREADS = int32(opts.NTHREADS);
        else                                NTHREADS = int32(0);                        end
    end
end

% Run cover tree findwithin
if kNN > size(X,2), kNN = size(X,2); end

if isa(X,'single')
    nn = findnearestMEX(CoverTree,X,Xq,int32(kNN),NTHREADS);
elseif isa(X,'double')
    nn = findnearestMEXD(CoverTree,X,Xq,int32(kNN),NTHREADS);
else
    warning('\n covertree_nnsearch: invalid class %s for X',class(X));
end

% Reorganize output
idxs    = nn.indices+1;
dists   = nn.distances;
opts.NTHREADS = NTHREADS;

return