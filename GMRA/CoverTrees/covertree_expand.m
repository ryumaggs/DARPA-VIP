function [expTree,exp_nIdxs] = covertree_expand ( Tree, Opts )

% MM_DEBUG: does not work well when many points at 0 distance from each other

if nargin<2, Opts = []; end

if ~isfield(Opts,'ExtRange')    Opts.ExtRange = 'min';  end
%if ~isfield(Opts,'MaxLevel')    Opts.MaxLevel = Inf;    end

% Set up the new tree
max_level           = max(Tree.levels(:,1));
changed_levels      = [];
if max_level == intmax('int32'),
    changed_levels = Tree.levels(:,1)==intmax('int32');
    Tree.levels(changed_levels,1) = max(Tree.levels(Tree.levels(:,1)<intmax('int32'),1))+1;
    Opts.ExtRange = 'min';
end
[~,root_level]      = covertree_get_root( Tree );
min_level           = min(Tree.levels(:,1));
max_level           = max(Tree.levels(:,1));

% Allocate space for new tree: need to compute its size
nIdxs = cell(max_level-min_level+1,1);
expTreelen = 0;
for j = min_level:max_level
    nIdxs{j-min_level+1}    = find(Tree.levels(:,1)==j)'-1;                                                                     % Get the nodes at scale j
    expTreelen              = expTreelen+(max_level-j+1)*length( nIdxs{j-min_level+1} );
end

% Put the root first
expTree             = Tree;
expTree.levels      = nan(expTreelen,5);
exp_nIdxs           = nan(expTreelen,1);
exp_nIdxs(1)        = 0;
% Go level by level and add the new vertices to the expanded tree
cur_Treeidx         = 1;
init_cur_expTreeIdx = size(Tree.levels,1)+1;
cur_expTreeIdx      = init_cur_expTreeIdx;

cur_children_idx = 1;
% Go through all the levels in the tree
for j = min_level:max_level
    j_rel                               = j-min_level+1;
    curupdate_Idx                       = nIdxs{j_rel}+1;                                       % Add to expTree the nodes at level j_rel of Tree
    exp_nIdxs(curupdate_Idx)            = nIdxs{j_rel};
    expTree.levels(curupdate_Idx,1)     = Tree.levels(nIdxs{j_rel}+1,1);
    if j==root_level,
        expTree.levels(curupdate_Idx,2) = -1;
    end
    
    % Go through each of the nodes at level j_rel. Copy each of them to the bottom level, enlarging the tree, and at the same time
    % rewiring the children of the node so that they have the deepest possible parent (which is a copy of the node under consideration some scales down)
    for k = 1:length(curupdate_Idx)
        % Get the children of the node
        children_nIdxs                       = covertree_get_children( Tree,exp_nIdxs(curupdate_Idx(k)) );
        %         % Remove children which are below the maximal scale
        %         if Opts.MaxLevel<Inf,
        %             children_nIdxs = children_nIdxs(Tree.levels(children_nIdxs+1,1)<=Opts.MaxLevel+root_level);
        %         end
        switch( Opts.ExtRange )
            case 'min'
                if ~isempty(children_nIdxs)
                    scales_to_add            = j:min([max(Tree.levels(children_nIdxs+1,1)),max_level]);
                else
                    scales_to_add            = j;
                end
            case 'max'
                scales_to_add                = j:max_level;
        end
        for p = 1:length(scales_to_add)                                        % Loop through the scales to be added from the current node downward
            if scales_to_add(p)==j,
                cur_addidx                   = curupdate_Idx(k);               % Update the current node at the current level
            elseif scales_to_add(p)==j+1,                                      % Update the "extended" node; its parent is an original node we had already places in expTree.
                cur_addidx                   = cur_expTreeIdx;
                expTree.levels(cur_addidx,2) = curupdate_Idx(k)-1;
                cur_expTreeIdx               = cur_expTreeIdx + 1;
            else
                cur_addidx                   = cur_expTreeIdx;                 % Update the "extended" node; its parent is an "extended" node
                expTree.levels(cur_addidx,2) = cur_expTreeIdx-2;               % It is -2 because -1 for the parent, and then -1 to convert to node index which is 0 based
                cur_expTreeIdx               = cur_expTreeIdx+1;
            end
            expTree.levels(cur_addidx,1)     = scales_to_add(p);
            children_to_be_linked            = find(Tree.levels(children_nIdxs+1,1)==scales_to_add(p)+1); % Find the children linked to the copy of the node at this scale
            expTree.levels(cur_addidx,4)     = cur_children_idx-1;              % Add these as children of the copy of the node at this scale
            if p<length(scales_to_add),
                expTree.levels(cur_addidx,3) = length(children_to_be_linked)+1; % If not at the bottom scale, add a new child which is a copy of the current node
            else
                expTree.levels(cur_addidx,3) = length(children_to_be_linked);
            end
            expTree.levels(children_nIdxs(children_to_be_linked)+1,2) = cur_addidx-1;
            if cur_addidx >= init_cur_expTreeIdx,
                % Debug only:  if ~isnan(exp_nIdxs(cur_addidx)) && exp_nIdxs(cur_addidx)~=exp_nIdxs(curupdate_Idx(k)), fprintf('\n\n WARNING: exp_nIdxs(%d) already has value %d, writing now %d.\n',cur_addidx,exp_nIdxs(cur_addidx),exp_nIdxs(curupdate_Idx(k))); end
                exp_nIdxs(cur_addidx)        = exp_nIdxs(curupdate_Idx(k));
            end
            if expTree.levels(cur_addidx,3)>0,
                if p<length(scales_to_add),                                    % Copy these children in the global children list in the new tree
                    expTree.levels(cur_children_idx:cur_children_idx+expTree.levels(cur_addidx,3)-1,5) = [children_nIdxs(children_to_be_linked);cur_expTreeIdx-1];
                else                                                           % Copy these children in the global children list in the new tree
                    expTree.levels(cur_children_idx:cur_children_idx+expTree.levels(cur_addidx,3)-1,5) = children_nIdxs(children_to_be_linked);
                end
                cur_children_idx = cur_children_idx+expTree.levels(cur_addidx,3);
            end
        end
    end
    cur_Treeidx = cur_Treeidx + length(curupdate_Idx);
end


%% Clean up
expTree.levels(cur_children_idx,5) = 0;                                        % Add the root note as the last entry in the children column

%% Actual size of tree level structure
actual_size       = find(isnan(exp_nIdxs), 1 )-1;
if ~isempty(actual_size),
    expTree.levels    = expTree.levels(1:actual_size,:);
    expTree.levels    = int32(expTree.levels);
    exp_nIdxs         = exp_nIdxs(1:actual_size);
    exp_nIdxs         = int32(exp_nIdxs);
end
%size(expTree.levels,1),

if ~isempty(changed_levels)
    expTree.levels(changed_levels,1) = intmax('int32');
end

return





