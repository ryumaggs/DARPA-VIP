function [Projections, tangentialCorrections] = IGWT2(gMRA, XGWT)

% Inverse Geometric Wavelet Transform
%
% Input:
%   gMRA        : structure of wavelet bases and translations
%   XGWT        : output from FGWT_combined
%
% Output:
%  Projections           : D-N-J matirx of recovered data at all scales
%  TangentialCorrections : D-N-J array of tangential corrections

%% Initialization
Projections             = [];
tangentialCorrections   = [];

if isempty(XGWT), fprintf('\n There is no IGWT to be computed.\n'); return; end

N                       = numel(XGWT.leafNodeLabels);
J                       = max(gMRA.Scales); % number of scales
nLeafNodes              = numel(gMRA.LeafNodes);

Projections             = zeros(gMRA.opts.AmbientDimension, N, J); % wavelets at all J scales

if gMRA.opts.addTangentialCorrections && nargout>1
    tangentialCorrections = zeros(gMRA.opts.AmbientDimension, N, J);
end

%% Pre-processing
if gMRA.opts.verbose > 0, fprintf('\n\t IGWT: pre-processing...'); end
Timing.preprocessing                                = cputime;

leaveswithPts                                       = unique(XGWT.leafNodeLabels);
%leaveswithPts_idxs                                  = ismember(gMRA.LeafNodes,leaveswithPts);
[sortedLeafNodeLabels,sortedLeafNodeLabels_idxs]    = sort(XGWT.leafNodeLabels);
sortedLeafNodeLabels_idxs                           = uint32(sortedLeafNodeLabels_idxs);
sortedLeafNodeLabels_counts                         = uint32(histc(sortedLeafNodeLabels,unique(sortedLeafNodeLabels)));
sortedLeafNodeLabels_cumcounts                      = [0;uint32(cumsum(sortedLeafNodeLabels_counts))];

Timing.preprocessing                                = cputime-Timing.preprocessing;
if gMRA.opts.verbose > 0, fprintf('done. (%.5f secs)',Timing.preprocessing); end

%% Go through the leafnodes and invert along the path to the root
if gMRA.opts.verbose > 0, fprintf('\n\t IGWT: computing inverse...'); end
Timing.inverting                                     = cputime;

for i = 1:length(leaveswithPts)
    net         = leaveswithPts(i);
    iRow        = find(gMRA.LeafNodes==net);
    netPts      = sortedLeafNodeLabels_idxs(sortedLeafNodeLabels_cumcounts(i)+1:sortedLeafNodeLabels_cumcounts(i+1));
    nPts        = sortedLeafNodeLabels_counts(i);
    
    if nPts == 0, continue; end
    
    cur_sum     = 0;
    x_tmp_prev  = 0;
    x_matj      = zeros(gMRA.opts.AmbientDimension,nPts,J);
    j_max       = gMRA.Scales(net);                                                                                             % Number of scales involved
    chain       = dpath(gMRA.cp, net);                                                                                          % Path to root
    
    for j = j_max:-1:1
        if ~isempty(gMRA.WavConsts{chain(j)})
            x_tmp = gMRA.WavConsts{chain(j)};
        else
            x_tmp = 0;
        end
        if ~isempty(gMRA.WavBases{chain(j)}) && ~isempty(XGWT.CelWavCoeffs{iRow,j})
            x_tmp = bsxfun(@plus,gMRA.WavBases{chain(j)}'*XGWT.CelWavCoeffs{iRow,j},x_tmp);
        end
        
        if gMRA.opts.addTangentialCorrections && j>1                                                                            % x_mat is corrected by adding tangential corrections
            if j<j_max,
                cur_sum = bsxfun(@plus,x_tmp_prev,cur_sum);
                x_TC = gMRA.ScalBasis{chain(j-1)}'*(gMRA.ScalBasis{chain(j-1)}*cur_sum);
                x_tmp = bsxfun(@minus,x_tmp,x_TC);
                if nargout>1,      
                    if ~isempty(x_TC),
                        if size(x_TC,2)==length(netPts)
                            tangentialCorrections(:,netPts, j) = x_TC;
                        else
                            tangentialCorrections(:,netPts, j) = repmat(x_TC,[1,length(netPts)]);
                        end
                    end
                end
            end
        end
        if size(x_tmp,2)==size(x_matj,2),
            x_matj(:,:,j)   = x_tmp;
        end
        x_tmp_prev      = x_tmp;
    end
    
    Projections(:,netPts,:) = cumsum(x_matj,3);
end

Timing.inverting    = cputime-Timing.inverting;
if gMRA.opts.verbose > 0, fprintf('done. (%.5f secs)',Timing.inverting); end

return