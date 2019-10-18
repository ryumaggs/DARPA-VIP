function scales = compute_scales(cp)

% Computes the scales of the nodes in a tree encoded by cp. The root node has scale 1.

nAllNets = numel(cp);
currentNodes = 1:nAllNets;

scales = zeros(nAllNets,1);
while any(currentNodes)    
    nonzeroNodes = find(currentNodes > 0);
    scales(nonzeroNodes) = scales(nonzeroNodes) + 1;
    currentNodes(nonzeroNodes) = cp(currentNodes(nonzeroNodes));
end

return