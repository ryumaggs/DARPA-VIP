function GMRA_VisualizePartition( X, gMRA, Partitions )

if ~iscell(Partitions), lPartitions{1} = Partitions; Partitions = lPartitions; end

for k = 1:length(Partitions)
    figure
    labels = zeros(size(X,2),1);
    for l = 1:length(Partitions{k})
        labels(gMRA.PointsInNet{Partitions{k}(l)})= l;
    end
    scatter3( X(1,:),X(2,:),X(3,:),20,labels,'filled' );
    title(sprintf('Partition %d, with %d components',k,length(Partitions{k})));
end

return