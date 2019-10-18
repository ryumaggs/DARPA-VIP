function node_idxs = covertree_get_leaves( Tree )

% Nodes with no children are leaves, provided that they are at level at most (level_root+covertree_opts.numlevel), 
% because below that level the covertree has not been constructed

levels     = Tree.levels(:,1);
%level_root = min(levels);
max_level  = max(levels);

node_idxs = find( ((Tree.levels(:,3)==0) & (levels<max_level)) | (levels>=max_level) );

% node_idxs = find( ((Tree.levels(:,3)==0) & (levels<=level_root+Tree.outparams(3))) ...
%                 | (levels>=level_root+Tree.outparams(3)) );
            
node_idxs = node_idxs-1;
            
return