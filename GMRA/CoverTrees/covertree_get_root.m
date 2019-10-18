function [idx,level] = covertree_get_root( Tree )

% Returns the (1-based) index of the root of the covertree, and its level

idx     = find(Tree.levels(:,2)==-1,1);                                                                                         % There should be only one root
level   = Tree.levels(idx,1);

return;
