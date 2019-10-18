function [means,stds,mins,maxs] = covertree_getleveldiffs( CoverTree )

curidx = 1;

means = zeros(size(CoverTree.levels,1),1);
stds = zeros(size(CoverTree.levels,1),1);
mins = zeros(size(CoverTree.levels,1),1);
maxs = zeros(size(CoverTree.levels,1),1);

for k = 1:size(CoverTree.levels,1)
    difflevels = single(CoverTree.levels(covertree_get_children( CoverTree, k-1 )+1,1 ))-single(CoverTree.levels(k,1));
    if ~isempty(difflevels)
        means(curidx)   = mean(difflevels);
        stds(curidx)    = std(difflevels);
        mins(curidx)    = min(difflevels);
        maxs(curidx)    = max(difflevels);
        curidx          = curidx + 1;
    end
end

means(curidx:end) = [];
stds(curidx:end) = [];
mins(curidx:end) = [];
maxs(curidx:end) = [];

return