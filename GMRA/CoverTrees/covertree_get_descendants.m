function descendant_list = covertree_get_descendants( Tree, node, allowed_children)

% Terrible code as far as performance goes, TBD

descendant_list = zeros(size(Tree.levels,1),1);
descendant_count = 0;
tbd = node;

while ~isempty(tbd),
    children_list       = covertree_get_children(Tree,tbd(1));
    if nargin>=3,
        children_list   = intersect(children_list,allowed_children);
    end;
    n_children_list     = length(children_list);
    descendant_list(descendant_count:descendant_count+children_list-1) = children_list;
    descendant_count    = descendant_count + n_children_list;
    tbd                 = [tbd(2:end);children_list];                                   % Add the children to the list of nodes tbd
end;

descendant_list = descendant_list(1:descendant_count);

% if nargin>=3,
%     % Keep only the descendants at scale j_0 or larger (i.e. finer)
%     descendant_list = descendant_list( Tree.levels(descendant_list,1)>=j_0 );
% end;

return;