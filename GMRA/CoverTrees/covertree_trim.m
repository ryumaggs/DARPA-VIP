function [cp,cp_idxs,pc,CoverTree_new] = covertree_trim( CoverTree, Opts )

%
% function [cp,cp_idxs,pc,CoverTree_new] = covertree_trim( CoverTree, Opts )
%
% IN:
%   CoverTree : covertree data structure as returned from covertree (or covertree_build)
%   [Opts]    : structure of options:
%       [TrimType]  : value in {'Scale','Size'}. Decides what variables to use to decide whether to trim a node or not
%       [TrimValue] : threshold on TrimType property of a node to decide whether to trim or not
%
% OUT:
%   cp              : tree structure flattened to vector of parents
%   cp_idxs         : cp_idxs(i) is the index of the covertree node in CoverTree corresponding to the i-th node in the new tree
%   pc              : #pts times 3 matrix encoding the children of each node
%   CoverTree_new   : new covertree
%
%

% (c) Mauro Maggioni, Duke University


if nargin<2,    Opts = struct('TrimType','Scale','TrimValue',max(CoverTree.levels(:,1))-1 );      end

switch( Opts.TrimType )
    case 'Scale'
        cp_idxs = find(CoverTree.levels(:,1)<min(CoverTree.levels(:,1))+Opts.TrimValue);
        [~,cp]  = ismember(CoverTree.levels(cp_idxs,2)+1,cp_idxs);
    case 'Size'
        sizes   = covertree_sizedescendants( CoverTree );
        cp_idxs = find(sizes>=Opts.TrimValue);
        [~,cp]  = ismember(CoverTree.levels(cp_idxs,2)+1,cp_idxs);
    otherwise
        warning('\n WARNING:covertree_trim: unknown TrimType %s',Opts.TrimType);
        cp      = CoverTree.levels(:,2)+1;
        cp_idxs = 1:size(CoverTree.levels,1);
end

if nargout>=3                                                                                                % Compute the parent->children map if requested
    [~,invmap] = ismember(1:size(CoverTree.levels,1),cp_idxs);
    cp_idxs_binary(size(CoverTree.levels,1),1) = false;
    cp_idxs_binary(cp_idxs) = true;
    pc                      = zeros(length(cp),3);
    if nargout>=4
        CoverTree_new       = zeros(length(cp),5);
    end
    cur_idx                 = 1;
    CoverTree_new = CoverTree;
    CoverTree_new.levels = zeros(length(cp),5,'int32');
    for k = 1:length(cp)
        pc(k,1)                             = cur_idx;
        children_tmp                        = covertree_get_children(CoverTree,cp_idxs(k)-1)+1;                 % Get the children (and change their index to 1-based)
        children_cp_idxs                    = invmap(children_tmp(cp_idxs_binary(children_tmp)));
        pc(k,2)                             = length(children_cp_idxs);
        pc(cur_idx:(cur_idx+pc(k,2)-1),3)   = children_cp_idxs;
        if nargout>=4
            CoverTree_new.levels(k,1)                           = CoverTree.levels(cp_idxs(k),1);
            CoverTree_new.levels(k,2)                           = cp(k)-1;
            CoverTree_new.levels(k,3)                           = pc(k,2);
            CoverTree_new.levels(k,4)                           = pc(k,1)-1;
            CoverTree_new.levels(cur_idx:(cur_idx+pc(k,2)-1),5) = children_cp_idxs-1;
        end
        cur_idx                             = cur_idx+pc(k,2);
    end
end

return
