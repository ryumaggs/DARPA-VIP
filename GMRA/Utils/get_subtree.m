function gMRA = get_subtree(gMRA, flags)

%
gMRA.flags = flags;

%%
keptNodes = find(flags==1);

inverseMap = zeros(1, length(gMRA.cp));
inverseMap(keptNodes) = 1:length(keptNodes);

new_cp = gMRA.cp(keptNodes');
innerNodes = find(new_cp>0);
new_cp(innerNodes) = inverseMap(new_cp(innerNodes));

gMRA.cp = new_cp;
gMRA.IniLabels = inverseMap(gMRA.IniLabels);
gMRA.Scales = compute_scales(gMRA.cp);
gMRA.isaleaf = leafnodes(gMRA.cp)';
gMRA.LeafNodes = find(gMRA.isaleaf);

gMRA.PointsInNet = gMRA.PointsInNet(keptNodes);
gMRA.Radii = gMRA.Radii(keptNodes);
gMRA.Sizes = gMRA.Sizes(keptNodes);
gMRA.Centers = gMRA.Centers(keptNodes);
gMRA.ScalFuns = gMRA.ScalFuns(keptNodes);
gMRA.Sigmas = gMRA.Sigmas(keptNodes);
gMRA.WavBases = gMRA.WavBases(keptNodes);
gMRA.WavConsts = gMRA.WavConsts(keptNodes);
gMRA.WavSingVals = gMRA.WavSingVals(keptNodes);

%
if gMRA.opts.sparsifying && ~gMRA.opts.orthogonalizing
    gMRA.Projections = gMRA.Projections(keptNodes);
end

%
if gMRA.opts.pruning,
    gMRA.epsEncodingCosts = gMRA.epsEncodingCosts(keptNodes);
end

%
if gMRA.opts.splitting
    gMRA.opts.mergePsiCapIntoPhi = false; 
    gMRA.WavDimsPsiCap = zeros(1,length(gMRA.cp));
    nonleafNodes = find(~gMRA.isaleaf);
    for i = 1:length(nonleafNodes)
        gMRA = splitting_WaveletBases(gMRA, nonleafNodes(i));
    end
end
