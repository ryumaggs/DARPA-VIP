function [D,idxs,dists,NNInfo] = distancefcn_covertree( pt_to, pt_from, covertree_dest, distancefcn, opts )

%
% function [D,idxs,dists,NNInfo] = distancefcn_covertree( pt_to, pt_from, covertree_dest, distancefcn, opts )
%
% Compute the kNN nearest neighbors of the points pt_from in the set pt_to, using a
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
%   D               : full matrix N_1xN_2 distances
%   idxs            : cell array of indices into pt_to, of length N_2
%   dists           : cell array of distances to pt_to, of length N_2
%   NNInfo          : structure return by nrsearch
%
% 
% (c) Mauro Maggioni
%

% Test molecular distances
%   X_1 = [1;0;0;2;0;0]; X_2 = 1/sqrt(2)*[0;1;0;1;0;0]; [D,idxs,dists,NNInfo] = distancefcn_covertree( X_2, X_1, [], 7, struct('data_classname',2,'ReturnD', true) );
%



if nargin<5 || isempty(opts),               opts = struct();                           end
if nargin<2,                                pt_from = pt_to;                           end
if nargin>=3 && ~isempty(covertree_dest),   opts.NNInfo.CoverTree = covertree_dest;    end
if nargin>=4 && ~isempty(distancefcn),      opts.distancefcn = distancefcn;            end
if ~isfield(opts,'kNN') || opts.kNN==0,     opts.kNN = size(pt_to,2);                  end
if ~isfield(opts,'data_classname'),         opts.data_classname = int32(0);            end
if ~isfield(opts,'ReturnD'),                opts.ReturnD = true;                       end

opts.ReturnAsArrays = false;
opts.distancefcn    = int32(opts.distancefcn);
opts.CoverTreeOpts  = struct('theta',0.5,'numlevels',int32(10),'minlevel',int32(0),'BLOCKSIZE',int32(1024),'distancefcn',int32(opts.distancefcn),'classname',int32(opts.data_classname));
if isfield(opts,'NTHREADS'), opts.CoverTreeOpts.NTHREADS = opts.NTHREADS;              end

[~,idxs,dists,NNInfo] = nrsearch( pt_to, pt_from, opts.kNN, 0, opts );

if isempty(pt_from),
    pt_from = pt_to;
end

for i = 1:length(idxs)
    goodidxs = find(idxs{i}>0);
    idxs{i}  = idxs{i}(goodidxs);
    dists{i} = dists{i}(goodidxs);
end

if opts.ReturnD
    D = zeros(size(pt_to,2),size(pt_from,2));
    
    for i = 1:size(D,2)
        D(idxs{i},i) = dists{i};
    end
else
    D = [];
end

return