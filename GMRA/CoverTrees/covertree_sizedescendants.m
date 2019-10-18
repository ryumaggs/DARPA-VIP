function sizes = covertree_sizedescendants( Tree )

% function sizes = covertree_sizedescendants( CoverTree )
%
% IN:
%   CoverTree : a covertree data struture
%
% OUT:
%   sizes     : a column vector where the k-th entry is the 1+number of descendants of the k-th node in CoverTree
%

% Initialize variables
sizes = ones(size(Tree.levels,1),1,'uint32');

% Get the maximum level in the tree
levels   = Tree.levels(:,1);
maxLevel = max(levels);
if maxLevel == intmax('int32'),
    maxLevel = max(Tree.levels(Tree.levels(:,1)<intmax('int32'),1));;
end
minLevel = min(levels);

% Go from the bottom of the tree up, updating the size of the parent of each node at each level
for k = maxLevel:-1:minLevel                                                                                                    % Go through levels
    idxs = find(levels==k);
    for i = 1:length(idxs)                                                                                                      % Go through nodes at current level
        parent_nidx = Tree.levels(idxs(i),2)+1;                                                                            % Get parent
        if parent_nidx>0                                                                                                        % Update count of parent (if current node is not root)
            sizes(parent_nidx) = sizes(parent_nidx)+sizes(idxs(i)); 
        end
    end
end

return