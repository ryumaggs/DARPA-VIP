function [DistMat,idxs,dists,NNInfo] = distancefcn_KNN_covertree( pt_to, pt_from, covertree_dest, distancefcn, opts )

%
% function [D,idxs,dists,NNInfo] = distancefcn_KNN_covertree( pt_to, pt_from, covertree_dest, distancefcn, opts )
%
% IN:
%   pt_to           : DxN_1 matrix of points to which distances need to be computed
%   [pt_from]       : DxN_2 matrix of points from which distances need to be computed. Default: pt_to.
%   [distancefcn]   : distance function to use. Default: 0 (Euclidean).
%   [covetree_dest] : covertree for destination points. Default: will be constructed
%   [opts]          : structure of options:
%                       [kNN]   : number of nearest neighbors to compute. Default: N_1
%                       [data_classname] : data type. 0 for vectors, 1 for images, 2 for molecular states. Default: 0 (vectors).
%                       [ReturnD]        : return matrix of distances. Default: true.
%
% OUT:
%   DistMat         : sparse N_1xN_2 matrix of distances between each point in pt_from and its kNN's in pt_to. Only computed if opts.ReturnD
%   idxs            : opts.kNN x N_2 matrix of indices into pt_to, for each point in pt_from
%   dists           : opts.kNN x N_2 matrix of distances to pt_to, for each point in pt_from
%   NNInfo          : structure with options, covertree, etc.. TBD
%
% Test with distancefcn_KNN_covertree_test
%
% (c) Mauro Maggioni
%



if nargin<5 || isempty(opts),               opts = struct();                            end
if nargin<2 || isempty(pt_from),            pt_from = pt_to;                            end
if nargin>=3 && ~isempty(covertree_dest),   NNInfo.CoverTree = covertree_dest;
                                       else NNInfo.CoverTree = [];                      end
if nargin>=4 && ~isempty(distancefcn),      opts.distancefcn = distancefcn;             end
if ~isfield(opts,'kNN') || opts.kNN==0,     opts.kNN = size(pt_to,2);                   end
if ~isfield(opts,'data_classname'),         opts.data_classname = int32(0);             end
if ~isfield(opts,'ReturnD'),                opts.ReturnD = false;                       end

opts.ReturnAsArrays = false;
opts.distancefcn    = int32(opts.distancefcn);
opts.CoverTreeOpts  = struct('theta',0.5,'numlevels',int32(10),'minlevel',int32(0),'BLOCKSIZE',int32(1024),'distancefcn',int32(opts.distancefcn),'classname',int32(opts.data_classname));
if isfield(opts,'NTHREADS'), opts.CoverTreeOpts.NTHREADS = opts.NTHREADS;              end
NNInfo.opts         = opts;

if isempty(NNInfo.CoverTree)
    tic
    [NNInfo.CoverTree,NNInfo.CoverTreeOpts] = covertree_build( pt_to,opts.CoverTreeOpts );                                      % Construct cover tree
    NNInfo.Timings.CoverTree = toc;
end
tic;
[idxs,dists,NNInfo.KNNOpts]                 = covertree_nnsearch( pt_to, NNInfo.CoverTree, pt_from, opts.kNN, opts.CoverTreeOpts );     % Find nearest neighbors with cover trees
NNInfo.Timings.covertree_nnsearch = toc;

if isempty(pt_from),
    pt_from = pt_to;
end

if opts.ReturnD                 %% TBD: change this to sparse format
    DistMat = sparsekNNDist2Mat( idxs,dists, size(pt_to,2), size(pt_from,2) );
else
    DistMat = [];
end

return