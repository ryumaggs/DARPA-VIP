function PartitionsAsLabels = GetPartitionsAsLabels( gMRA )

J = max(gMRA.Scales);

PartitionsAsLabels = zeros(length(gMRA.IniLabels),J);

leaves_idxs = get_partition_at_scale( gMRA, J );

for k = 1:length(leaves_idxs)
    cur_node = leaves_idxs(k);
    cur_idxs = gMRA.IniLabels==cur_node;    
    while cur_node>0
        PartitionsAsLabels(cur_idxs,gMRA.Scales(cur_node)) = cur_node;        
        cur_node = gMRA.cp(cur_node);
    end
end

return