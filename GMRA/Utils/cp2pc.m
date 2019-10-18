function pc = cp2pc( cp )

% Convert children to parent map to parent to children map
pc = zeros(length(cp),3);
cur_pcidx = 1;
for k = 1:length(cp)                                                                                        % Done the slow way
    idxs                                = find(cp==k);
    pc(k,1)                             = cur_pcidx;
    pc(k,2)                             = length(idxs);
    pc(cur_pcidx:cur_pcidx+pc(k,2)-1,3) = idxs;
    cur_pcidx                           = cur_pcidx+pc(k,2);
end

return