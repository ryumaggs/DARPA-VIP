function n_idxs = covertree_get_pathtoroot( Tree, node )

n_idxs = node;

while Tree.levels(node+1,2)~=-1,
    node = Tree.levels(node+1,2);
    n_idxs = [n_idxs,node];
end;

return;