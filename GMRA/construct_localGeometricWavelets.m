function nodeWavInfo = construct_localGeometricWavelets(X,gMRA,node,children)

% OUT:
%   nodeWavInfo structure with the following fields:
%       WavBases    : wavelet basis (->gMRA.WavBases{children(c)}
%       WavSingVals : wavelet singular values (->WavSingVals{children(c)}
%       WavConsts   : wavelet centers (->WavConsts{children(c)}

if nargin<3; children = find(gMRA.cp==node); end;

nChildren = numel(children);

%%
allFineBases = [vertcat(gMRA.ScalBasis{children})];

if ~gMRA.opts.orthogonalizing
    parentScalBasis = gMRA.ScalBasis{node};
else
    parentScalBasis = [gMRA.WavBases{[node get_ancestors(gMRA.cp, node)]}];                                                     % orthogonal matrix by induction
end

% if ~gMRA.opts.sparsifying                                                                                                     % TBD: this does not really make sense if Predictor is true
%     allFineBasesPerp = allFineBases - (allFineBases*parentScalBasis')* parentScalBasis;
% else
%     allFineBasesPerp = allFineBases - parentScalBasis * (parentScalBasis\allFineBases);                                       % TBD
% end

%% Compute the wavelet subspaces W_{j+1,k'}
wavDims = zeros(1,nChildren);
%col     = 0;
for c = 1:nChildren
    switch gMRA.opts.Predictor.PredictorType
        case 'none'
            %   Y2                               = allFineBasesPerp(col+1:col+size(gMRA.ScalBasis{children(c)}, 1),:);                     % Difference vectors children-parent projections
            if ~gMRA.opts.sparsifying                                                                                                       % TBD: this does not really make sense if Predictor is true
                Y                               = gMRA.ScalBasis{children(c)}-(gMRA.ScalBasis{children(c)}*parentScalBasis')*parentScalBasis;
            else
                Y                               = gMRA.ScalBasis{children(c)}-parentScalBasis*(parentScalBasis\gMRA.ScalBasis{children(c)});    % TBD?
            end
            %    norm(Y-Y2),
            [~,S,V]                         = randPCA(Y,min([min(size(Y)),gMRA.opts.MaxDim]));                                          % Use randomized PCA. Otherwise: [U,S] = svd(Y, 'econ');
            wavDims(c)                      = sum(diag(S)>gMRA.opts.threshold1);%*max(size(Y)) * eps(norm(S,inf)));                     % Wavelet dimension determined by threshold1
            if wavDims(c)>0
                nodeWavInfo.WavBases{c}     = V(:,1:wavDims(c));
                nodeWavInfo.WavSingVals{c}  = diag(S(1:wavDims(c), 1:wavDims(c)));
            else
                nodeWavInfo.WavBases{c}     = [];
                nodeWavInfo.WavSingVals{c}  = [];
            end            
            nodeWavInfo.WavConsts{c}        = gMRA.Centers{children(c)} - gMRA.Centers{node};
            if gMRA.opts.sparsifying
                nodeWavInfo.WavConsts{c}    = nodeWavInfo.WavConsts{c} - parentScalBasis*(parentScalBasis\nodeWavInfo.WavConsts{c});
            else
                nodeWavInfo.WavConsts{c}    = nodeWavInfo.WavConsts{c} - parentScalBasis'*(parentScalBasis*nodeWavInfo.WavConsts{c});
            end
            %% Splitting of the wavelet subspaces for the children nodes, if requested
            if gMRA.opts.splitting
                gMRA = splitting_WaveletBases(gMRA, node, children);
            end
            
            %% Sparsifying basis in the wavelet subspaces, if requested
            if gMRA.opts.sparsifying
                gMRA = sparsifying_WaveletBases(gMRA, node, children);                                                                      % TBD
            end
        %    col                             = col+size(gMRA.ScalBasis{children(c)}, 1);
        case 'orthogonal'
            Xparent     = parentScalBasis'*parentScalBasis*bsxfun(@minus,X(:,gMRA.PointsInNet{children(c)}),gMRA.Centers{node});
            Xchildren   = gMRA.ScalBasis{children(c)}'*gMRA.ScalBasis{children(c)}*bsxfun(@minus,X(:,gMRA.PointsInNet{children(c)}),gMRA.Centers{children(c)});
            [nodeWavInfo.Odiscr(children(c)),XparentO, nodeWavInfo.O{children(c)}] = procrustes(Xchildren',Xparent','Scaling',false);         % Transform Xparent to Xchildren
            Y           = Xchildren-XparentO';
            nodeWavInfo.WavSingVals{c}  = norm(Y,'fro')^2;
            nodeWavInfo.WavConsts{c}    = [];
            nodeWavInfo.WavBases{c}     = [];            
    end
end

return


nodeWavInfo2=nodeWavInfo;



% 
% 
% 
% 
% 
% function nodeWavInfo = construct_localGeometricWavelets(gMRA,node,children)
% 
% % OUT: 
% %   nodeWavInfo structure with the following fields:
% %       WavBases    : wavelet basis (->gMRA.WavBases{children(c)}
% %       WavSingVals : wavelet singular values (->WavSingVals{children(c)}
% %       WavConsts   : wavelet centers (->WavConsts{children(c)}
% 
% if nargin<3; children = find(gMRA.cp==node); end;
% 
nChildren = numel(children);

%%
allFineBases = [vertcat(gMRA.ScalBasis{children})];
% allFineBases = [];
% for c = 1:nChildren
%     allFineBases = [allFineBases gW.ScalBasis{children(c)}*diag(gW.Sigmas{children(c)}(1:size(gW.ScalBasis{children(c)},2)))];
% end

if ~gMRA.opts.orthogonalizing
    parentScalBasis = gMRA.ScalBasis{node};  
else   
    parentScalBasis = [gMRA.WavBases{[node get_ancestors(gMRA.cp, node)]}]; % orthogonal matrix by induction
end
    
if ~gMRA.opts.sparsifying
    allFineBasesPerp = allFineBases - (allFineBases*parentScalBasis')* parentScalBasis;
else
    allFineBasesPerp = allFineBases - parentScalBasis * (parentScalBasis\allFineBases);
end

%% Compute the wavelet subspaces W_{j+1,k'}
wavDims = zeros(1,nChildren);
col     = 0;
for c = 1:nChildren    
    Y                               = allFineBasesPerp(col+1:col+size(gMRA.ScalBasis{children(c)}, 1),:);                       % Difference vectors children-parent projections
    [~,S,V]                         = randPCA(Y,min([min(size(Y)),gMRA.opts.MaxDim]));                                          % Use randomized PCA. Otherwise: [U,S] = svd(Y, 'econ');
    wavDims(c)                      = sum(diag(S)>gMRA.opts.threshold1);%*max(size(Y)) * eps(norm(S,inf)));                     % Wavelet dimension determined by threshold1
    if wavDims(c)>0
        nodeWavInfo.WavBases{c}     = V(:,1:wavDims(c));   
        nodeWavInfo.WavSingVals{c}  = diag(S(1:wavDims(c), 1:wavDims(c)));
    else
        nodeWavInfo.WavBases{c}     = [];
        nodeWavInfo.WavSingVals{c}  = [];
    end
    
    nodeWavInfo.WavConsts{c}        = gMRA.Centers{children(c)} - gMRA.Centers{node};
    if gMRA.opts.sparsifying
        nodeWavInfo.WavConsts{c}    = nodeWavInfo.WavConsts{c} - parentScalBasis*(parentScalBasis\nodeWavInfo.WavConsts{c});
    else
        nodeWavInfo.WavConsts{c}    = nodeWavInfo.WavConsts{c} - parentScalBasis'*(parentScalBasis*nodeWavInfo.WavConsts{c});
    end        
    col                             = col+size(gMRA.ScalBasis{children(c)}, 1);    
end

%% Splitting of the wavelet subspaces for the children nodes, if requested
if gMRA.opts.splitting
    gMRA = splitting_WaveletBases(gMRA, node, children);
end

%% Sparsifying basis in the wavelet subspaces, if requested
if gMRA.opts.sparsifying
    gMRA = sparsifying_WaveletBases(gMRA, node, children);                                                                      % TBD
end

for c = 1:length(nodeWavInfo.WavBases)
if norm(nodeWavInfo.WavBases{c}-nodeWavInfo2.WavBases{c})>0 || ...
    norm(nodeWavInfo.WavConsts{c}-nodeWavInfo2.WavConsts{c})>0
    keyboard
end
end


return
