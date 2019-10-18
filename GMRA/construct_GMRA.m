function gMRA = construct_GMRA( X, gMRA )

% 
% function gMRA = construct_GMRA( gMRA )
%
% Constructs a GMRA structure
% 

%% Construct the leaf nodes first, at the finest scale
if gMRA.opts.verbose>0,         fprintf('\n\t..constructing scaling spaces for leaf nodes...'); end
gMRA    = construct_vNets_leafNodes( X, gMRA );

%% Go bottom-up, from fine to coarser scales
J     = max(gMRA.Scales);
for j = J-1:-1:1   
    nodes           = find((gMRA.Scales==j) & (~gMRA.isaleaf));                                                                 % Find nodes at scale j
    cp              = gMRA.cp;                                                                                                  % Get the parent
    PointsInNet     = gMRA.PointsInNet;                                                                                         % Get the points in the node
    NewPointsInNet  = cell(1,length(nodes));
    nodeInfo        = cell(1,length(nodes));
    children        = cell(1,length(nodes));
    if gMRA.opts.verbose>0,         fprintf('\n\t..constructing scale %d...',j); end
    if gMRA.opts.parallel
        for k  = 1:length(nodes)
            children{k}                 = find(cp==nodes(k));                                                                   % Go through the nodes at scale j
            NewPointsInNet{k}           = [PointsInNet{children{k}}];                                                           % Get the children
        end
        parfor k = 1:length(nodes)
            nodeInfo{k}                 = local_SVD_analysis(X,gMRA,nodes(k),NewPointsInNet{k});                                  % Construct the scaling functions at node (j,k)
        end
    else        
        for k = 1:length(nodes)                                                                                                 % Go through the nodes at scale j
            children{k}                 = gMRA.pc(gMRA.pc(nodes(k),1):(gMRA.pc(nodes(k),1)+gMRA.pc(nodes(k),2)-1),3);
            NewPointsInNet{k}           = [PointsInNet{children{k}}];                                                           % Update points in node
            nodeInfo{k}                 = local_SVD_analysis(X,gMRA,nodes(k),NewPointsInNet{k});                                % Construct the scaling functions at node (j,k)
        end
    end
    for k = 1:length(nodes)                                                                                                     % Copy the results from local SVD analysis to GMRA structure
        gMRA.PointsInNet{nodes(k)}      = NewPointsInNet{k};
        gMRA.Sizes(nodes(k))            = nodeInfo{k}.Size;
        gMRA.Centers{nodes(k)}          = nodeInfo{k}.Center;
        gMRA.Radii(nodes(k))            = nodeInfo{k}.Radius;
        gMRA.Sigmas{nodes(k)}           = nodeInfo{k}.Sigmas;
        gMRA.ScalBasis{nodes(k)}        = nodeInfo{k}.ScalBasis';
        if isfield(nodeInfo{k},'epsEncodingCost')
            gMRA.epsEncodingCost(nodes(k))  = nodeInfo{k}.epsEncodingCost;
        end
        if isfield(nodeInfo{k},'Projections')
            gMRA.Projections{nodes(k)}      = nodeInfo{k}.Projections;
        end
    end    
    if gMRA.opts.ComputeWavelets,
        if gMRA.opts.parallel
            parfor k = 1:length(nodes)                                                                                          % Go through the nodes at scale j
                nodeWavInfo{k}              = construct_localGeometricWavelets(X, gMRA,nodes(k),children{k});                      % Construct the wavelets (j+1,children(j,k))
            end
        else
            for k = 1:length(nodes)                                                                                             % Go through the nodes at scale j
                nodeWavInfo{k}              = construct_localGeometricWavelets(X, gMRA,nodes(k),children{k});                      % Construct the wavelets (j+1,children(j,k))
            end
        end
        for k = 1:length(nodes)                                                                                                 % Copy the geometric wavelets just constructed to GMRA structure
            for i = 1:length(children{k})
                gMRA.WavBases{children{k}(i)}       = nodeWavInfo{k}.WavBases{i}';
                gMRA.WavSingVals{children{k}(i)}    = nodeWavInfo{k}.WavSingVals{i};
                gMRA.WavConsts{children{k}(i)}      = nodeWavInfo{k}.WavConsts{i};
                if isfield(nodeWavInfo{k},'O')
                    gMRA.O{children{k}(i)}          = nodeWavInfo{k}.O;
                    gMRA.Odiscr{children{k}(i)}     = nodeWavInfo{k}.Odiscr;
                else
                    gMRA.O{children{k}(i)}          = [];
                end
            end
        end
    else
        for k = 1:length(nodes)                                                                                                 % Copy the geometric wavelets just constructed to GMRA structure
            for i = 1:length(children{k})
                gMRA.WavBases{children{k}(i)}      = [];
                gMRA.WavSingVals{children{k}(i)}   = [];
                gMRA.WavConsts{children{k}(i)}     = [];
            end
        end        
    end
end


%% Compute other statistics for the tree
tmp_Radii_AvePerScale   = zeros(1,J);
tmp_Radii               = gMRA.Radii;
tmp_Scales              = gMRA.Scales;
for j = 1:J
    tmp_Radii_AvePerScale(j) = mean(tmp_Radii(tmp_Scales==j));
end
gMRA.Radii_AvePerScale  = tmp_Radii_AvePerScale;

% Wavelet bases at the root are just the scaling functions
roots                   = find(gMRA.cp==0);
gMRA.WavBases(roots)    = gMRA.ScalBasis(roots);
gMRA.WavConsts(roots)   = gMRA.Centers(roots);
gMRA.WavSingVals(roots) = gMRA.Sigmas(roots);

return
