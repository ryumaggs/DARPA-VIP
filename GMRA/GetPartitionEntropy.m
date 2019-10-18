function [PartitionEntropy , TreeEntropy] = GetPartitionEntropy( gMRA, XGWT, Partition )

%% function entropy = GetPartitionEntropy( gMRA, XGWT, Partition )
%% 
%% compute the entropy of a thresholded tree

PartitionEntropy = sum(gMRA.Radii(Partition.nodes).^2);

nodes   = [];
for i = 1:size(XGWT.Cel_cpidx,1)
    nodes   = [nodes XGWT.Cel_cpidx(i,1:Partition.leafScale(i))];
end

nodes       = unique(nodes(nodes>0));
TreeEntropy = sum(gMRA.Radii(nodes).^2);



% entropy = 0;
% for i = 1:size(XGWT.Cel_cpidx,1)
%     uniquev = XGWT.Cel_cpidx(i,1:gMRA.Scales(Partition(i)));
%     uniquev = uniquev(uniquev>0);
%     entropy = entropy + sum(gMRA.Radii(uniquev).^2);
% end

return
