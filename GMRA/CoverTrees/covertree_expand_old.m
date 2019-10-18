function [expTree,exp_nIdxs] = covertree_expand ( Tree )

% Set up the new tree
[~,root_level]      = covertree_get_root( Tree );
max_level           = max(Tree.levels(:,1)-root_level);

% Allocate space for new tree: need to compute its size
nIdxs = cell(max_level-root_level+1,1);
expTreelen = 0;
for j = root_level:max_level
    % Count number of nodes at scale j
    nIdxs{j-root_level+1} = find(Tree.levels(:,1)==j)'-1;
    expTreelen = 2*expTreelen+length( nIdxs{j-root_level+1} );
end;


% Put the root first
expTree             = Tree;
expTree.levels      = -1*ones(expTreelen,5,'int32');
exp_nIdxs           = nan(expTreelen,1);

%expTree.levels(1,1:2) = Tree.levels(1,1:2);
exp_nIdxs(1)         = 0;
% Go level by level and add the new vertices to the expanded tree
cur_Treeidx         = 1;
init_cur_expTreeIdx = size(Tree.levels,1)+1;
cur_expTreeIdx      = init_cur_expTreeIdx;
prev_levelidxs{1}   = [];

cur_children_idx = 1;
% Go through all the scales in the tree
for j = root_level:max_level
    j_rel                                                       = j-root_level+1;
    % Add to expTree the nodes at level j_rel of Tree
    curupdate_Idx                                               = cur_Treeidx:cur_Treeidx+length(nIdxs{j_rel})-1;
    exp_nIdxs(curupdate_Idx)                                    = nIdxs{j_rel};
    old_nodes                                                   = find(nIdxs{j_rel}<init_cur_expTreeIdx);
    expTree.levels(exp_nIdxs(curupdate_Idx(old_nodes))+1,1:2)   = Tree.levels(nIdxs{j_rel}(old_nodes)+1,1:2);
    
    % Now look at the children, see if the current node needs to be copied one level down, and rewire children as needed
    for k = 1:length(curupdate_Idx)
        % Get the children of the node
        children_nIdxs       = covertree_get_children( Tree,exp_nIdxs(curupdate_Idx(k)) );
        % Copy the node down the finest scale
        scales_to_add = j:max_level;
        for p = 1:length(scales_to_add)
            if scales_to_add(p)==j,
                % We are going to update the current node at the current level
                cur_addidx                   = exp_nIdxs(curupdate_Idx(k))+1;
            elseif scales_to_add(p)==j+1,
                % We are going to update the "extended" node and its parent is an original node we had already places in expTree.
                cur_addidx                   = cur_expTreeIdx;
                expTree.levels(cur_addidx,2) = exp_nIdxs(curupdate_Idx(k));
                cur_expTreeIdx               = cur_expTreeIdx + 1;
            else
                % We are going to update the "extended" node and its parent is an "extended" node
                cur_addidx                   = cur_expTreeIdx;
                expTree.levels(cur_addidx,2) = cur_expTreeIdx-2;                % It is -2 because -1 for the parent, and then -1 to convert to node index which is 0 based
                cur_expTreeIdx               = cur_expTreeIdx+1;
            end
            expTree.levels(cur_addidx,1)     = scales_to_add(p);
            % Find the children that are to be linked to the copy of the node at this scale
            children_to_be_linked = find(Tree.levels(children_nIdxs+1,1)==scales_to_add(p)+1);
            % Add these as children of the copy of the node at this scale
            expTree.levels(cur_addidx,4)     = cur_children_idx-1;
            % If not at the bottom scale, add a child which is a copy of the current node
            if p<length(scales_to_add),
                expTree.levels(cur_addidx,3)     = length(children_to_be_linked)+1;
            else
                expTree.levels(cur_addidx,3)     = length(children_to_be_linked);
            end;
            if cur_addidx >= init_cur_expTreeIdx,
                if exp_nIdxs(cur_addidx) ~= -1 && exp_nIdxs(cur_addidx)~=exp_nIdxs(curupdate_Idx(k)),
                    warning(sprintf('exp_nIdxs(%d) already has value %d, writing now %d.',cur_addidx,exp_nIdxs(cur_addidx),exp_nIdxs(curupdate_Idx(k))));
                end;
                exp_nIdxs(cur_addidx)        = exp_nIdxs(curupdate_Idx(k));
            end;
            if expTree.levels(cur_addidx,3)>0,
                if p<length(scales_to_add),
                    % Copy these children in the global children list in the new tree
                    expTree.levels(cur_children_idx:cur_children_idx+expTree.levels(cur_addidx,3)-1,5) = [children_nIdxs(children_to_be_linked);cur_expTreeIdx-1];
                    exp_nIdxs(cur_expTreeIdx)   = exp_nIdxs(curupdate_Idx(k));
                else
                    % Copy these children in the global children list in the new tree
                    expTree.levels(cur_children_idx:cur_children_idx+expTree.levels(cur_addidx,3)-1,5) = children_nIdxs(children_to_be_linked);
                end;
                cur_children_idx = cur_children_idx+expTree.levels(cur_addidx,3);
            end;
        end;
    end;
    cur_Treeidx = cur_Treeidx + length(curupdate_Idx);
end;


%% Clean up
% Add the root note as the last entry in the children column
expTree.levels(cur_children_idx,5) = 1;

% Actual size of tree level structure
actual_size       = find(isnan(exp_nIdxs), 1 )-1;
expTree.levels    = expTree.levels(1:actual_size,:);
exp_nIdxs         = exp_nIdxs(1:actual_size);
exp_nIdxs         = int32(exp_nIdxs);

return









%%

% % %
% % %
% % %
% % % %         % Get the children of every copied node
% % % %         if expIdxs(cur_updateidx(k)) < init_cur_expTreeidx,
% % % %             children_idxs       = covertree_get_children( Tree,expIdxs(cur_updateidx(k)) );
% % % %         else
% % % %             children_idxs       = covertree_get_children( expTree,expIdxs(cur_updateidx(k)) );
% % % %         end
% % % %         % Children that are at most one level down from current node, stay
% % % %         children_that_stay  = find(Tree.levels(children_idxs,1)<=j+1);
% % % %         if ~isempty(children_that_stay),
% % % %             expTree.levels(cur_children_idx:cur_children_idx+length(children_that_stay)-1,5) = children_idxs(children_that_stay);
% % % %             expTree.levels(expIdxs(cur_updateidx(k)),3) = cur_children_idx;
% % % %             expTree.levels(expIdxs(cur_updateidx(k)),4) = length(children_that_stay);
% % % %             cur_children_idx=cur_children_idx+length(children_that_stay);
% % % %         else
% % % %             expTree.levels(expIdxs(cur_updateidx(k)),3) = 0;
% % % %             expTree.levels(expIdxs(cur_updateidx(k)),4) = 0;
% % % %         end;
% % % %         % Children that are more than one level down from current node, require copy of current node one scale down, and rewiring
% % % %         children_that_dont_stay = find(Tree.levels(children_idxs,1)>j+1);
% % % %         % expParent_idx = cur_expTreeidx(find(expIdxs(cur_expTreeidx)==expIdxs(cur_updateidx(k))));
% % % %         if ~isempty(children_that_dont_stay),
% % % %             % Add children that do not stay to global children list
% % % %             expTree.levels(cur_children_idx:cur_children_idx+length(children_that_dont_stay)-1,5) = children_idxs(children_that_dont_stay);
% % % %             % If the node has not been copied one level down, copy it one level down
% % % %             expParent_idx = find(expIdxs(init_cur_expTreeidx:cur_expTreeidx)==expIdxs(cur_updateidx(k)));
% % % %             if ~isempty(expParent_idx),
% % % %                 expParent_idx   = init_cur_expTreeidx+expParent_idx-1;
% % % %             else
% % % %                 expParent_idx   = cur_expTreeidx;
% % % %                 % Add added node to list of nodes at the next scale
% % % %                 idxs{j_rel+1} = union(idxs{j_rel+1},expParent_idx);%expIdxs(cur_updateidx(k)));
% % % %                 cur_expTreeidx  = cur_expTreeidx + 1;
% % % %             end;
% % % %             expTree.levels(expParent_idx,1) = j+1;
% % % %             expTree.levels(expParent_idx,2) = expIdxs(cur_updateidx(k));
% % % %             % Assign the children to this new node
% % % %             expTree.levels(expParent_idx,3) = cur_children_idx;
% % % %             expTree.levels(expParent_idx,4) = length(children_that_dont_stay);
% % % %             expIdxs(expParent_idx)          = expIdxs(cur_updateidx(k));
% % % %
% % % %             % Update the current position in the list of current children and list of nodes added to tree
% % % %             cur_children_idx                 = cur_children_idx+length(children_that_dont_stay);
% % % %         end;
% % % end;
% % %
% % %
% % %
% % % %     if j_rel>1,
% % % %         prev_levelidxs{j_rel}                                                   = (cur_Treeidx-length(idxs{j_rel})):(cur_Treeidx+length(prev_levelidxs{j_rel-1})-1);
% % % %         cur_Treeidx                                                          = cur_Treeidx + length(prev_levelidxs{j_rel-1});
% % % %     else
% % % %         prev_levelidxs{j_rel}                                                   = cur_updateidx;
% % % %         cur_Treeidx                                                          = cur_Treeidx + length(cur_updateidx);
% % % %     end;
% % % cur_Treeidx = cur_Treeidx + length(cur_updateidx);
% % % end;
% % %
% % %
% % % %% Clean up
% % % % Add the root note as the last entry in the children column
% % % expTree.levels(cur_children_idx,5) = 1;
% % %
% % % % Actual size of tree level structure
% % % actual_size = min(find(isnan(expIdxs)))-1;
% % % expTree.levels = expTree.levels(1:actual_size,:);
% % % expIdxs = expIdxs(1:actual_size);
% % %
% % % return
% % %
% % %
% % %
% % %
% % %
% % %
% % %
% % %
% % %
% % %


































% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % % function [expTree,expIdxs] = covertree_expand ( Tree )
% % % %
% % % % % Set up the new tree
% % % % [root,root_level] = covertree_get_root( Tree );
% % % % max_level = max(Tree.levels(:,1)-root_level);
% % % %
% % % % % Allocate space for new tree: need to compute its size
% % % % idxs = cell(max_level-root_level+1,1);
% % % % expTreelen = 0;
% % % % for j = root_level:max_level
% % % %     % Count number of nodes at scale j
% % % %     idxs{j-root_level+1} = find(Tree.levels(:,1)==j)';
% % % %     expTreelen = 2*expTreelen+length( idxs{j-root_level+1} );
% % % % end;
% % % %
% % % %
% % % % % Put the root first
% % % % expTree             = Tree;
% % % % expTree.levels      = nan(expTreelen,5);
% % % % expIdxs             = nan(size(expTree.levels,1),1);
% % % %
% % % % %expTree.levels(1,1:2) = Tree.levels(1,1:2);
% % % % expIdxs(1)     = 1;
% % % % % Go level by level and add the new vertices to the expanded tree
% % % % cur_Treeidx      = 1;
% % % % cur_expTreeidx   = sum([idxs{:}])+1;
% % % % prev_levelidxs{1}   = [];
% % % %
% % % % cur_children_idx = 1;
% % % % for j = root_level:max_level
% % % %     j_rel                                                                   = j-root_level+1;
% % % %     % Copy over the nodes at level j of Tree
% % % %     cur_updateidx                                                           = cur_Treeidx:cur_Treeidx+length(idxs{j_rel})-1;
% % % %     expTree.levels(cur_updateidx,1:2)                                       = Tree.levels(idxs{j_rel},1:2);
% % % %     expIdxs(cur_updateidx)                                                  = idxs{j_rel};
% % % %     if j_rel>1,
% % % %         % Now add nodes from the previous level of the expTree
% % % %         cur_Treeidx                                                          = cur_Treeidx + length(idxs{j_rel});
% % % %         cur_updateidx                                                       = cur_Treeidx:cur_Treeidx+length(prev_levelidxs{j_rel-1})-1;
% % % %         expTree.levels(cur_updateidx,1:2)                                   = expTree.levels(prev_levelidxs{j_rel-1},1:2);
% % % %         % Update parents: these are the same points as those at the previous level
% % % %         expTree.levels(cur_updateidx,2)                                     = prev_levelidxs{j_rel-1};
% % % %         % The points "copied" from the previous scale now live at the new scale...
% % % %         expTree.levels(cur_updateidx,1)                                         = j;
% % % %         % ...and they have no children
% % % %         expTree.levels(cur_updateidx,3)                                         = 0;
% % % %         % Save the absolute index (uid) of these points
% % % %         expIdxs(cur_updateidx)                                                  = expIdxs(prev_levelidxs{j_rel-1});
% % % %     else
% % % %         cur_updateidx                                                       = cur_Treeidx;
% % % %     end;
% % % %     % Now rewire the children
% % % %     for k = 1:length(cur_updateidx)
% % % %         % Get the children of every copied node
% % % %         children_idxs = covertree_get_children( Tree,expIdxs(cur_updateidx(k)) );
% % % %         children_that_stay = find(Tree.levels(children_idxs,1)<=j);
% % % %         if ~isempty(children_that_stay),
% % % %             expTree.levels(cur_children_idx:cur_children_idx+length(children_that_stay)-1,5) = children_idxs(children_that_stay);
% % % %             expTree.levels(expIdxs(cur_updateidx(k)),3) = cur_children_idx;
% % % %             expTree.levels(expIdxs(cur_updateidx(k)),4) = length(children_that_stay);
% % % %             cur_children_idx=cur_children_idx+length(children_that_stay);
% % % %         else
% % % %             expTree.levels(expIdxs(cur_updateidx(k)),3) = 0;
% % % %             expTree.levels(expIdxs(cur_updateidx(k)),4) = 0;
% % % %         end;
% % % %         children_that_dont_stay = find(Tree.levels(children_idxs,1)>j+1);
% % % %         if ~isempty(children_that_dont_stay),
% % % %             expTree.levels(cur_children_idx:cur_children_idx+length(children_that_dont_stay)-1,5) = children_idxs(children_that_dont_stay);
% % % %             expTree.levels(cur_Treeidx,3) = cur_children_idx;
% % % %             expTree.levels(cur_Treeidx,4) = length(children_that_dont_stay);
% % % %             cur_children_idx                 = cur_children_idx+length(children_that_dont_stay);
% % % %         else
% % % %             expTree.levels(cur_Treeidx,3) = 0;
% % % %             expTree.levels(cur_Treeidx,4) = 0;
% % % %         end;
% % % %     end;
% % % %
% % % %     if j_rel>1,
% % % %         prev_levelidxs{j_rel}                                                   = (cur_Treeidx-length(idxs{j_rel})):(cur_Treeidx+length(prev_levelidxs{j_rel-1})-1);
% % % %         cur_Treeidx                                                          = cur_Treeidx + length(prev_levelidxs{j_rel-1});
% % % %     else
% % % %         prev_levelidxs{j_rel}                                                   = cur_updateidx;
% % % %         cur_Treeidx                                                          = cur_Treeidx + length(cur_updateidx);
% % % %     end;
% % % % end;
% % % %
% % % %
% % % % return
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %     function [expTree,expIdxs] = covertree_expand ( Tree )
% % % %
% % % %         % Set up the new tree
% % % %         [root,root_level] = covertree_get_root( Tree );
% % % %         max_level = max(Tree.levels(:,1)-root_level);
% % % %
% % % %         % Allocate space for new tree: need to compute its size
% % % %         idxs = cell(max_level-root_level+1,1);
% % % %         expTreelen = 0;
% % % %         for j = root_level:max_level
% % % %             % Count number of nodes at scale j
% % % %             idxs{j-root_level+1} = find(Tree.levels(:,1)==j)';
% % % %             expTreelen = 2*expTreelen+length( idxs{j-root_level+1} );
% % % %         end;
% % % %
% % % %
% % % %         % Put the root first
% % % %         expTree             = Tree;
% % % %         expTree.levels      = nan(expTreelen,5);
% % % %         expIdxs             = nan(expTreelen,1);
% % % %
% % % %         %expTree.levels(1,1:2) = Tree.levels(1,1:2);
% % % %         expIdxs(1)          = 1;
% % % %         % Go level by level and add the new vertices to the expanded tree
% % % %         cur_Treeidx         = 1;
% % % %         init_cur_expTreeidx = size(Tree.levels,1)+1;
% % % %         cur_expTreeidx      = init_cur_expTreeidx;
% % % %         prev_levelidxs{1}   = [];
% % % %
% % % %         %expTree.levels(1:size(Tree.levels,1),:) = Tree.levels;
% % % %         cur_children_idx = 1;
% % % %         for j = root_level:max_level
% % % %             j_rel                                                                   = j-root_level+1;
% % % %             % Add the nodes at level j of Tree to expTree
% % % %             cur_updateidx                                           = cur_Treeidx:cur_Treeidx+length(idxs{j_rel})-1;
% % % %             expIdxs(cur_updateidx)                                  = idxs{j_rel};
% % % %             old_nodes                                               = find(idxs{j_rel}<init_cur_expTreeidx);
% % % %             expTree.levels(expIdxs(cur_updateidx(old_nodes)),1:2)   = Tree.levels(idxs{j_rel}(old_nodes),1:2);
% % % %
% % % %
% % % %
% % % %
% % % %             %
% % % %             %     if j_rel>1,
% % % %             %         % Copy nodes from the previous level of the expTree
% % % %             %         %        cur_Treeidx                                                          = cur_Treeidx + length(idxs{j_rel});
% % % %             %         cur_expTreeidx                                                       = cur_expTreeidx:cur_expTreeidx+length(prev_levelidxs{j_rel-1})-1;
% % % %             %         expTree.levels(cur_expTreeidx,1:2)                                   = expTree.levels(prev_levelidxs{j_rel-1},1:2);
% % % %             %         % Update parents: these are the same points as those at the previous level
% % % %             %         expTree.levels(cur_expTreeidx,2)                                     = prev_levelidxs{j_rel-1};
% % % %             %         % The points "copied" from the previous scale now live at the new scale...
% % % %             %         expTree.levels(cur_expTreeidx,1)                                         = j;
% % % %             %         % ...and they have no children
% % % %             %         expTree.levels(cur_expTreeidx,3)                                         = 0;
% % % %             %         % Save the absolute index (uid) of these points
% % % %             %         expIdxs(cur_expTreeidx)                                                  = expIdxs(prev_levelidxs{j_rel-1});
% % % %             %     end;
% % % %             % Now look at the children, see if the current node needs to be copied one level down, and rewire children as needed
% % % %             for k = 1:length(cur_updateidx)
% % % %                 % Get the children of every copied node
% % % %                 if expIdxs(cur_updateidx(k)) < init_cur_expTreeidx,
% % % %                     children_idxs       = covertree_get_children( Tree,expIdxs(cur_updateidx(k)) );
% % % %                 else
% % % %                     children_idxs       = covertree_get_children( expTree,expIdxs(cur_updateidx(k)) );
% % % %                 end
% % % %                 % Children that are at most one level down from current node, stay
% % % %                 children_that_stay  = find(Tree.levels(children_idxs,1)<=j+1);
% % % %                 if ~isempty(children_that_stay),
% % % %                     expTree.levels(cur_children_idx:cur_children_idx+length(children_that_stay)-1,5) = children_idxs(children_that_stay);
% % % %                     expTree.levels(expIdxs(cur_updateidx(k)),3) = cur_children_idx;
% % % %                     expTree.levels(expIdxs(cur_updateidx(k)),4) = length(children_that_stay);
% % % %                     cur_children_idx=cur_children_idx+length(children_that_stay);
% % % %                 else
% % % %                     expTree.levels(expIdxs(cur_updateidx(k)),3) = 0;
% % % %                     expTree.levels(expIdxs(cur_updateidx(k)),4) = 0;
% % % %                 end;
% % % %                 % Children that are more than one level down from current node, require copy of current node one scale down, and rewiring
% % % %                 children_that_dont_stay = find(Tree.levels(children_idxs,1)>j+1);
% % % %                 % expParent_idx = cur_expTreeidx(find(expIdxs(cur_expTreeidx)==expIdxs(cur_updateidx(k))));
% % % %                 if ~isempty(children_that_dont_stay),
% % % %                     % Add children that do not stay to global children list
% % % %                     expTree.levels(cur_children_idx:cur_children_idx+length(children_that_dont_stay)-1,5) = children_idxs(children_that_dont_stay);
% % % %                     % If the node has not been copied one level down, copy it one level down
% % % %                     expParent_idx = find(expIdxs(init_cur_expTreeidx:cur_expTreeidx)==expIdxs(cur_updateidx(k)));
% % % %                     if ~isempty(expParent_idx),
% % % %                         expParent_idx   = init_cur_expTreeidx+expParent_idx-1;
% % % %                     else
% % % %                         expParent_idx   = cur_expTreeidx;
% % % %                         % Add added node to list of nodes at the next scale
% % % %                         idxs{j_rel+1} = union(idxs{j_rel+1},expParent_idx);%expIdxs(cur_updateidx(k)));
% % % %                         cur_expTreeidx  = cur_expTreeidx + 1;
% % % %                     end;
% % % %                     expTree.levels(expParent_idx,1) = j+1;
% % % %                     expTree.levels(expParent_idx,2) = expIdxs(cur_updateidx(k));
% % % %                     % Assign the children to this new node
% % % %                     expTree.levels(expParent_idx,3) = cur_children_idx;
% % % %                     expTree.levels(expParent_idx,4) = length(children_that_dont_stay);
% % % %                     expIdxs(expParent_idx)          = expIdxs(cur_updateidx(k));
% % % %
% % % %                     % Update the current position in the list of current children and list of nodes added to tree
% % % %                     cur_children_idx                 = cur_children_idx+length(children_that_dont_stay);
% % % %                 end;
% % % %             end;
% % % %
% % % %
% % % %
% % % %             %     if j_rel>1,
% % % %             %         prev_levelidxs{j_rel}                                                   = (cur_Treeidx-length(idxs{j_rel})):(cur_Treeidx+length(prev_levelidxs{j_rel-1})-1);
% % % %             %         cur_Treeidx                                                          = cur_Treeidx + length(prev_levelidxs{j_rel-1});
% % % %             %     else
% % % %             %         prev_levelidxs{j_rel}                                                   = cur_updateidx;
% % % %             %         cur_Treeidx                                                          = cur_Treeidx + length(cur_updateidx);
% % % %             %     end;
% % % %             cur_Treeidx = cur_Treeidx + length(cur_updateidx);
% % % %         end;
% % % %
% % % %
% % % %         return
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %
% % % %         function [expTree,expIdxs] = covertree_expand ( Tree )
% % % %
% % % %             % Set up the new tree
% % % %             [root,root_level] = covertree_get_root( Tree );
% % % %             max_level = max(Tree.levels(:,1)-root_level);
% % % %
% % % %             % Allocate space for new tree: need to compute its size
% % % %             idxs = cell(max_level-root_level+1,1);
% % % %             expTreelen = 0;
% % % %             for j = root_level:max_level
% % % %                 % Count number of nodes at scale j
% % % %                 idxs{j-root_level+1} = find(Tree.levels(:,1)==j)';
% % % %                 expTreelen = 2*expTreelen+length( idxs{j-root_level+1} );
% % % %             end;
% % % %
% % % %
% % % %             % Put the root first
% % % %             expTree             = Tree;
% % % %             expTree.levels      = nan(expTreelen,5);
% % % %             expIdxs             = nan(size(expTree.levels,1),1);
% % % %
% % % %             %expTree.levels(1,1:2) = Tree.levels(1,1:2);
% % % %             expIdxs(1)     = 1;
% % % %             % Go level by level and add the new vertices to the expanded tree
% % % %             cur_Treeidx      = 1;
% % % %             cur_expTreeidx   = sum([idxs{:}])+1;
% % % %             prev_levelidxs{1}   = [];
% % % %
% % % %             cur_children_idx = 1;
% % % %             for j = root_level:max_level
% % % %                 j_rel                                                                   = j-root_level+1;
% % % %                 % Copy over the nodes at level j of Tree
% % % %                 cur_updateidx                                                           = cur_Treeidx:cur_Treeidx+length(idxs{j_rel})-1;
% % % %                 expTree.levels(cur_updateidx,1:2)                                       = Tree.levels(idxs{j_rel},1:2);
% % % %                 expIdxs(cur_updateidx)                                                  = idxs{j_rel};
% % % %                 if j_rel>1,
% % % %                     % Now add nodes from the previous level of the expTree
% % % %                     cur_Treeidx                                                          = cur_Treeidx + length(idxs{j_rel});
% % % %                     cur_updateidx                                                       = cur_Treeidx:cur_Treeidx+length(prev_levelidxs{j_rel-1})-1;
% % % %                     expTree.levels(cur_updateidx,1:2)                                   = expTree.levels(prev_levelidxs{j_rel-1},1:2);
% % % %                     % Update parents: these are the same points as those at the previous level
% % % %                     expTree.levels(cur_updateidx,2)                                     = prev_levelidxs{j_rel-1};
% % % %                     % The points "copied" from the previous scale now live at the new scale...
% % % %                     expTree.levels(cur_updateidx,1)                                         = j;
% % % %                     % ...and they have no children
% % % %                     expTree.levels(cur_updateidx,3)                                         = 0;
% % % %                     % Save the absolute index (uid) of these points
% % % %                     expIdxs(cur_updateidx)                                                  = expIdxs(prev_levelidxs{j_rel-1});
% % % %                 else
% % % %                     cur_updateidx                                                       = cur_Treeidx;
% % % %                 end;
% % % %                 % Now rewire the children
% % % %                 for k = 1:length(cur_updateidx)
% % % %                     % Get the children of every copied node
% % % %                     children_idxs = covertree_get_children( Tree,expIdxs(cur_updateidx(k)) );
% % % %                     children_that_stay = find(Tree.levels(children_idxs,1)<=j);
% % % %                     if ~isempty(children_that_stay),
% % % %                         expTree.levels(cur_children_idx:cur_children_idx+length(children_that_stay)-1,5) = children_idxs(children_that_stay);
% % % %                         expTree.levels(expIdxs(cur_updateidx(k)),3) = cur_children_idx;
% % % %                         expTree.levels(expIdxs(cur_updateidx(k)),4) = length(children_that_stay);
% % % %                         cur_children_idx=cur_children_idx+length(children_that_stay);
% % % %                     else
% % % %                         expTree.levels(expIdxs(cur_updateidx(k)),3) = 0;
% % % %                         expTree.levels(expIdxs(cur_updateidx(k)),4) = 0;
% % % %                     end;
% % % %                     children_that_dont_stay = find(Tree.levels(children_idxs,1)>j+1);
% % % %                     if ~isempty(children_that_dont_stay),
% % % %                         expTree.levels(cur_children_idx:cur_children_idx+length(children_that_dont_stay)-1,5) = children_idxs(children_that_dont_stay);
% % % %                         expTree.levels(cur_Treeidx,3) = cur_children_idx;
% % % %                         expTree.levels(cur_Treeidx,4) = length(children_that_dont_stay);
% % % %                         cur_children_idx                 = cur_children_idx+length(children_that_dont_stay);
% % % %                     else
% % % %                         expTree.levels(cur_Treeidx,3) = 0;
% % % %                         expTree.levels(cur_Treeidx,4) = 0;
% % % %                     end;
% % % %                 end;
% % % %
% % % %                 if j_rel>1,
% % % %                     prev_levelidxs{j_rel}                                                   = (cur_Treeidx-length(idxs{j_rel})):(cur_Treeidx+length(prev_levelidxs{j_rel-1})-1);
% % % %                     cur_Treeidx                                                          = cur_Treeidx + length(prev_levelidxs{j_rel-1});
% % % %                 else
% % % %                     prev_levelidxs{j_rel}                                                   = cur_updateidx;
% % % %                     cur_Treeidx                                                          = cur_Treeidx + length(cur_updateidx);
% % % %                 end;
% % % %             end;
% % % %
% % % %
% % % %             return