function covertree_buildmultiscalepartition( X, Tree )

C = 1/8;

% For every scale j, find ball of radius C2^{-j} around each point
j_min = min( Tree.levels(:,1) );
j_max = max( Tree.levels(:,2) );
N     = size(Tree.levels,1);

Cjkcores = cell(N,1);

for j = 1:(j_max-j_min+1)
    r    = Tree.radii(j);
    idxs{j} = find( Tree.levels(:,1)==j_min+j-1 );
    Cjkcores{idxs{j}} = covertree_rangesearch( X, Tree, X(:,idxs{j}), C*r );
end

% Create a tree connecting the cores
cur_idx = 1;
cur_children_idx = 1;
for j = 1:(j_max-j_min+1)
    for k = 1:length(idxs{j})
        levels(idxs{j}(k),1) = j;
        levels(idxs{j}(k),2) = 0;
        levels(idxs{j}(k),[3,4]) = 
    end
end

return