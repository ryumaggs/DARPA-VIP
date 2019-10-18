function children_n_idxs = covertree_get_children( Tree, node )

% Get the children of a given node in the covertree Tree, node is the (0-based) index of the node
% The children_n_idxs returned are the (0-based) indices of the children of node.

n_children      = Tree.levels(node+1,3);                                              % Find how many children the node has
n_startidx      = Tree.levels(node+1,4);                                              % Find where the list of children starts
children_n_idxs = Tree.levels((n_startidx+1):(n_startidx+n_children),5);              % These are the (0-based) indices of the children

return;