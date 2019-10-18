function DistMat = sparsekNNDist2Mat( idxs, dists, N_to, N_from )

count  = numel(idxs);
i_idxs = zeros(count,1);
j_idxs = zeros(count,1);

curidx = 1;
step   = size(idxs,1);
for j = 1:size(idxs,2)
    i_idxs(curidx:curidx+step-1) = idxs(:,j);
    j_idxs(curidx:curidx+step-1) = j;
    curidx = curidx + step;
end

DistMat = sparse(i_idxs,j_idxs,double(dists(:)),N_to,N_from,2*count);

return