function gMRA = construct_vNets_leafNodes( X,gMRA )

nLeafnodes      = length(gMRA.LeafNodes);
nodes           = zeros(1,nLeafnodes);
nodeInfo        = cell(1,nLeafnodes);
NewPointsInNet  = cell(1,nLeafnodes);
LeafNodes       = gMRA.LeafNodes;

if gMRA.opts.verbose > 0, fprintf('\n\t construct_vNets_leafNodes: pre-processing...'); end
Timing.preprocessing                                = cputime;

[sortedLeafNodeLabels,sortedLeafNodeLabels_idxs]    = sort(gMRA.IniLabels);
sortedLeafNodeLabels_idxs                           = uint32(sortedLeafNodeLabels_idxs);
sortedLeafNodeLabels_counts                         = uint32(histc(sortedLeafNodeLabels,unique(sortedLeafNodeLabels)));
sortedLeafNodeLabels_cumcounts                      = [0,uint32(cumsum(sortedLeafNodeLabels_counts))];

Timing.preprocessing                                = cputime-Timing.preprocessing;
if gMRA.opts.verbose > 0, fprintf('done. (%.5f secs)',Timing.preprocessing); end


if gMRA.opts.parallel
    parfor k = 1:length(sortedLeafNodeLabels_cumcounts)-1,
        nodes(k)                = LeafNodes(k);                                                                             % Get the index of the current node
        NewPointsInNet{k}       = sortedLeafNodeLabels_idxs(sortedLeafNodeLabels_cumcounts(k)+1:sortedLeafNodeLabels_cumcounts(k+1));
        nodeInfo{k}             = local_SVD_analysis(X,gMRA,nodes(k),NewPointsInNet{k});                                      % Perform local SVD analysis
    end
else
    for k = 1:length(sortedLeafNodeLabels_cumcounts)-1,
        nodes(k)                = LeafNodes(k);                                                                             % Get the index of the current node
        NewPointsInNet{k}       = sortedLeafNodeLabels_idxs(sortedLeafNodeLabels_cumcounts(k)+1:sortedLeafNodeLabels_cumcounts(k+1));
        nodeInfo{k}             = local_SVD_analysis(X,gMRA,nodes(k),NewPointsInNet{k});                                      % Perform local SVD analysis
    end
end
for k = 1:length(nodes)
    gMRA.PointsInNet{nodes(k)}  = NewPointsInNet{k};
    gMRA.Sizes(nodes(k))        = nodeInfo{k}.Size;
    gMRA.Centers{nodes(k)}      = nodeInfo{k}.Center;
    gMRA.Radii(nodes(k))        = nodeInfo{k}.Radius;
    gMRA.Sigmas{nodes(k)}       = nodeInfo{k}.Sigmas;
    gMRA.ScalBasis{nodes(k)}    = nodeInfo{k}.ScalBasis';
    
    if isfield(nodeInfo{k},'epsEncodingCost')   gMRA.epsEncodingCost(nodes(k))  = nodeInfo{k}.epsEncodingCost; end
    if isfield(nodeInfo{k},'Projections')       gMRA.Projections{nodes(k)}      = nodeInfo{k}.Projections;     end
end;

return;
