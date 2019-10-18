function [idxsflag,outerleaves] = GetSmallestSubTreeWithIdxs( cp, idxskeep, pc )

% function idxs = GetSmallestSubTreeWithIdxs( cp, idxskeep )
%
% IN:
%   cp          : vector encoding tree
%   idxskeep    : indices of the nodes of the tree to be kept, as a logical vector of length(cp)
%   [pc]        : parent->children structure associate to cp. Output of cp2pc(cp).
%
% OUT:
%   idxs        : indices of the nodes in the smallest subtree of cp containing all idxskeep==true nodes
%   outerleaves : indices of the nodes that are outer leaves
%
% (c) Mauro Maggioni

% MM: the implementation is not efficient, at least about a log(|cp|) off.

if nargin<3 || isempty(pc)
    pc = cp2pc( cp );
end

idxsflag    = idxskeep;

% for k = 1:length(idxskeep)                                                                                                      % Go through the nodes to be kept. Not the most efficient way of doing this, but alas!
%     if idxskeep(k)
%         pathtoroot                          = dpath(cp, k);
%         idxsflag( pathtoroot )              = true;                                                                             % Keep all the nodes in the path to the root from a node to be kept
%     end
% end
%
% idxsflag_orig = idxsflag;
% idxsflag = idxskeep;

idxs = find(idxskeep);

for k = 1:length(idxs)                                                                                                          % Go through the nodes to be kept, and fill up the tree from there.
    curnode = idxs(k);
    while true
        idxsflag(curnode) = true;
        curnode = cp(curnode);
        if curnode==0 || idxsflag(curnode)                                                                                      % Stop at root, of it a node already marked is met
            break
        end
    end
end

% if sum(idxsflag~=idxsflag_orig)
%     keyboard;
% end


if nargout>1
    outerleaves = zeros(size(idxskeep),'uint8');
    idxs = find(idxsflag);
    for k = 1:length(idxs)                                                                                                  % Pick the outer leaves
        children_idxs = pc(pc(idxs(k),1):(pc(idxs(k),1)+pc(idxs(k),2)-1),3);
        if isempty(children_idxs)
            outerleaves(idxs(k)) = true;
        else
            outerleaves(children_idxs(~idxsflag(children_idxs))) = true;
        end
    end
    outerleaves = find(outerleaves);
end

%% Check that we have a partition
if false
    leaves = find(pc(:,2)==0);
    for k = 1:length(leaves)
        lpath = dpath(cp,leaves(k));
        intersection = intersect(lpath,outerleaves);
        if length(intersection)~=1,
            warning('ops!')
            keyboard
        end
    end
end

return
