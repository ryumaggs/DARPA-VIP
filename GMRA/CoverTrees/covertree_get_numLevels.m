function numLevels = covertree_get_numLevels( Tree )

levels = Tree.levels(:,1);
levels(levels>=10000) = 0;

numLevels = range(levels)+1;

return;