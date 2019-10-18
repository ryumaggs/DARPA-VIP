function net = covertree_get_netatlevel( Tree, j )

net = find(Tree.levels(:,1)<=j);

return;