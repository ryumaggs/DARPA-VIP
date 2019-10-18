function idxs = get_partition_at_scale( GMRA, j )

%
% Returns a partition of the set using nodes at scale j, plus any additional leaf that may be needed at scale coarser than j.
% The return set of indices is indexing into the set of nodes of the GMRA tree.
% 

idxs = sort([find(GMRA.Scales == j);GMRA.LeafNodes(GMRA.Scales(GMRA.LeafNodes)<j)], 'ascend');

return