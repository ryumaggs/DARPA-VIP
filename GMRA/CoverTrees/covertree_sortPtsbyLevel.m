function sorted_idxs = covertree_sortPtsbyLevel( CoverTree )

[~,sorted_idxs] = sort(CoverTree.levels(:,1),'ascend');

return