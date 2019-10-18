function gMRA = construct_orthogonalGeometricWavelets( X, gMRA )

% build the tree, and compute svd only at the roots
root                = find(gMRA.cp==0);
gMRA                = construct_vNets_NoSVD(gMRA,root);       
nodeInfo            = local_SVD_analysis(X,gMRA,root,gMRA.PointsInNet{root});
gMRA.Sizes(root)    = nodeInfo.Size;
gMRA.Centers{root}  = nodeInfo.Center;
gMRA.Radii(root)    = nodeInfo.Radius;
gMRA.Sigmas{root}   = nodeInfo.Sigmas;
gMRA.ScalFuns{root} = nodeInfo.ScalFuns;
if isfield(nodeInfo,'epsEncodingCost')   gMRA.epsEncodingCost(root)  = nodeInfo.epsEncodingCost; end
if isfield(nodeInfo,'Projections')       gMRA.Projections{root}      = nodeInfo.Projections;     end

% initialization
gMRA.WavBases(root)     = gMRA.ScalFuns(root);
gMRA.WavConsts(root)    = gMRA.Centers(root);
gMRA.WavSingVals(root)  = gMRA.Sigmas(root);

% vector of indicators (whether to keep (1) or remove (-1) the node) 
flags = ones(1, numel(gMRA.cp)); 

%%
J = max(gMRA.Scales);
nAllNodes = length(gMRA.cp);

% scale 1
j = 1;
parentNodes = root; % set of nodes whose children will be processed (i.e., computing wavelet bases)

PointsInNet = gMRA.PointsInNet;

while j<J && ~isempty(parentNodes);        
    childrenNodes = zeros(1,nAllNodes); % collecting nodes for next round    
    for i = 1:length(parentNodes)        
        node = parentNodes(i);        
        if ~gMRA.isaleaf(node)            
            children = find(gMRA.cp==node);                        
            nodeInfo = cell(1,length(children));
            if gMRA.opts.parallel
                parfor c = 1:length(children)
                    nodeInfo{c} = local_SVD_analysis(X,gMRA,children(c),PointsInNet{children(c)});
                end
            else
                for c = 1:length(children)
                    nodeInfo{c} = local_SVD_analysis(X,gMRA,children(c),PointsInNet{children(c)});
                end
            end
            for c = 1:length(children)
                gMRA.Sizes(children(c))            = nodeInfo{c}.Size;
                gMRA.Centers{children(c)}          = nodeInfo{c}.Center;
                gMRA.Radii(children(c))            = nodeInfo{c}.Radius;
                gMRA.Sigmas{children(c)}           = nodeInfo{c}.Sigmas;
                gMRA.ScalFuns{children(c)}         = nodeInfo{c}.ScalFuns;
                if isfield(nodeInfo{c},'epsEncodingCost')   gMRA.epsEncodingCost(children(c))  = nodeInfo.epsEncodingCost; end
                if isfield(nodeInfo{c},'Projections')       gMRA.Projections{children(c)}      = nodeInfo.Projections;     end                                
            end
            nodeWavInfo = construct_localGeometricWavelets(gMRA,node,children);            
            for c = 1:length(children)
                gMRA.WavBases{children(c)}      = nodeWavInfo.WavBases{c};
                gMRA.WavSingVals{children(c)}   = nodeWavInfo.WavSingVals{c};
                gMRA.WavConsts{children(c)}     = nodeWavInfo.WavConsts{c};
            end;
            
            ancestors = get_ancestors(gMRA.cp, node);           
            for c = 1:length(children)                
                cumDict = [gMRA.WavBases{[children(c) node ancestors]}]; % orthogonal matrix by induction                
                Y = bsxfun(@minus, gMRA.X(:,gMRA.PointsInNet{children(c)}), gMRA.Centers{children(c)});                
                Y_proj = cumDict'*Y;
                approxErr = sum(Y.^2, 1) - sum((Y_proj).^2,1);                
                if strcmpi(gMRA.opts.errorType, 'absolute')
                    approxErr = sqrt(mean(approxErr));
                else
                    sigs = svd(Y_proj);
                    approxErr = sum(approxErr)/(sum(sigs.^2)+sum(approxErr));
                end                
                if approxErr <= gMRA.opts.precision % already within precision
                    offspring = get_offspring(gMRA.cp, children(c));
                    flags(offspring) = -1;
                    gMRA.IniLabels(gMRA.PointsInNet{children(c)}) = children(c);
                else
                    childrenNodes(children(c)) = 1;
                end                
            end            
        end        
    end    
    j = j+1;
    parentNodes = find(childrenNodes>0);    
end

%%
gMRA = get_subtree(gMRA, flags);

return