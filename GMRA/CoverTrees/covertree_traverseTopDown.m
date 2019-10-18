function [nodeList,parentList] = covertree_traverseTopDown( coverTree )

%
% function covertree_traverseTopDown( coverTree )
%
% IN:
%   coverTree   : a covertree
%
% OUT:
%   nodeList    : indices of nodes sorted from top (root) to bottom
%   parentList  : indices of the parent nodes of the nodes in nodeList (0 denotes no parent)
%
% (c) Mauro Maggioni, 2015

[~,nodeList]    = sort(coverTree.levels(:,1));
if nargout>1
    parentList      = coverTree.levels(nodeList,2)+1;
end

return