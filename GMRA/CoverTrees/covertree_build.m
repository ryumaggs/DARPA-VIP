function [CoverTree,opts] = covertree_build( X, opts )

%
% function [idxs,dists] = covertree_build( cX, opts )
%
% Wrapper for cover tree covertree construction function
%
% IN:
%   cX              : D by N matrix of N points in R^D
%   [opts]          : structure of options.
%                       [NTHREADS]      : number of threads for search. Default: number of cores.
%                       [distancfcn]    : index of distance function, int32.
%                       [classname]     : index of class name, int32.
%                       [numlevels]     : number of levels, int32.
%                       [minlevels]     : minimum level, int32.
%                       [theta]         : shrinkage factor from scale to scale, same class as X.
%                       [BLOCKSIZE]     : size of block of points per thread, int32.
%
% OUT:
%   CoverTree       : covertree data structure, with the following fields:
%                       theta   : dilation factor (0.5 is the default)
%                       radii   : row vector of radii at different levels
%                       levels  : n by 5 matrix. the i-th row contains [level,parent,#children,index first children,children idxs]
%                                 for the i-th point.
%                       outparams : 9 vector of parameters
%                       nclasstogetdist : number of distance computations
%
%
% (c) Copyright Duke University, 2014[
% Mauro Maggioni
% mauro@math.duke.edu
%
persistent NTHREADS

CoverTree = [];

MIN_N_FOR_MULTITHREADED  = 4096;

if nargin<2,                    opts = struct();                                        end
if isunix && ~ismac
    opts.NTHREADS = int32(0);
else
    if size(X,2)>MIN_N_FOR_MULTITHREADED
        if isfield(opts,'NTHREADS'),        opts.NTHREADS = int32(opts.NTHREADS);
        else                                opts.NTHREADS = int32(feature('numcores'));     end
    else
        if isfield(opts,'NTHREADS'),        opts.NTHREADS = int32(opts.NTHREADS);
        else                                opts.NTHREADS = int32(0);                       end
    end
end
if ~isfield(opts,'BLOCKSIZE'),      opts.BLOCKSIZE      = int32(2048);                  end
if ~isfield(opts,'distancefcn'),    opts.distancefcn    = int32(0);                     end
if ~isfield(opts,'classname'),      opts.classname      = int32(0);                     end
if ~isfield(opts,'theta'),          opts.theta          = single(0.5);                  end
if ~isfield(opts,'numlevels'),      opts.numlevels      = int32(100);                   end
if ~isfield(opts,'minlevel'),       opts.minlevel       = int32(0);                     end
if ~isfield(opts,'NTREES'),         opts.NTREES         = int32(1);                     end

if ~isempty(X)
    if isa(X,'single')
        opts.theta  = single(opts.theta);
        tic;
        CoverTree   = covertree( opts, X );
        Covertree.timing = toc;
    elseif isa(X,'double')
        opts.theta  = double(opts.theta);
        tic;
        CoverTree   = covertreeD( opts, X );
        Covertree.timing = toc;
    else
        warning('\n covertree_build: invalid class %s for X',class(X));
    end
end

return