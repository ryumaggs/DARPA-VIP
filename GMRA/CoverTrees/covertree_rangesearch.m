function [idxs,dists] = covertree_rangesearch( X, CoverTree, Xq, radius, opts )

%
% function [idxs,dists] = covertree_rangesearch( X, CoverTree, Xq, radius, opts )
%
% Wrapper for cover tree findwithinMEX function
%
% IN:
%   X               : D by N matrix of N points in R^D
%   CoverTree       : covertree on cX, constructed by covertree
%   Xq              : D by M matrix M query points in R^D
%   radius          : radius of search
%   [opts]          : structure of options.
%                       [NTHREADS]      : number of threads for search. Default: number of cores.
%
% OUT:
%   idxs            : cell array with the k-th cell containing the indices (into columns of cX) of the points within distance radius from cXq(:,k)
%   dists           : cell array with the distances
%   The above outputs are sorted in increasing distance.
%
%
% (c) Copyright Duke University, 2013
% Mauro Maggioni
% mauro@math.duke.edu
%

persistent NTHREADS

MIN_N_FOR_MULTITHREADED  = 1024;

if nargin<5,                    opts = struct();                            end

if isunix && ~ismac
    NTHREADS = int32(0);
else
    if size(X,2)>MIN_N_FOR_MULTITHREADED
        if isfield(opts,'NTHREADS')
            NTHREADS = int32(opts.NTHREADS);
        else
            NTHREADS = int32(feature('numcores'));
        end
    else
        NTHREADS = int32(0);
    end
end

if isa(X,'single')
    covertree.rangesearch = findwithinMEX(CoverTree,single(X),struct('distances',single(radius),'numlevels',int32(10)),single(Xq),NTHREADS);
elseif isa(X,'double')
    covertree.rangesearch = findwithinMEXD(CoverTree,double(X),struct('distances',double(radius),'numlevels',int32(10)),double(Xq),NTHREADS);
else
    warning('\n covertree_rangesearch: invalid class %s for X',class(X));
end

% Reorganize output
idxs    = cell(size(Xq,2),1);
dists   = cell(size(Xq,2),1);

for k = 1:length(idxs)
    idxs{k}     = covertree.rangesearch.indices(covertree.rangesearch.pi(k,2)+1:covertree.rangesearch.pi(k,2)+covertree.rangesearch.pi(k,1))+1;
    dists{k}    = covertree.rangesearch.distances(covertree.rangesearch.pi(k,2)+1:covertree.rangesearch.pi(k,2)+covertree.rangesearch.pi(k,1));
    goodidxs    = idxs{k}>0;
    idxs{k}     = idxs{k}(goodidxs);
    dists{k}    = dists{k}(goodidxs);
end;


return;
