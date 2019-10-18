function [ApproxError,Quantiles] = GWT_displayResults( X, gMRA, XGWT, imgOpts, Sigma, ApproxError, Quantiles)

[D, N] = size(X);

if nargin<4,                            imgOpts                 = struct();     end
if ~isfield(imgOpts, 'imageData'),      imgOpts.imageData       = false;        end
if ~isfield(imgOpts, 'isCompressed'),   imgOpts.isCompressed    = false;        end
if nargin<5 || isempty(Sigma),          Sigma                   = 0;            end

% Display the coefficents
GWT_DisplayCoeffs( gMRA, XGWT );

% Plot approximation error
if nargin<6
    ErrOpts = struct('norm',[2,inf],'relative',true,'quantiles',[0.25,0.5,0.75]);                                               % Default option for approximation error calculation
    [ApproxError,Quantiles] = GetApproximationErrors(X, gMRA, XGWT, ErrOpts );
end

DisplayApproximationErrors( gMRA,ApproxError.Absolute,Quantiles.Absolute,Sigma );
DisplayApproximationErrors( gMRA,ApproxError.Relative,Quantiles.Relative,Sigma );

% If image data, display some more stuff
if imgOpts.imageData,
    
    if isfield(imgOpts, 'sampleImage')
        i = imgOpts.sampleImage;
    else
        i = randsample(N,1);
    end
    
    leafNode = gMRA.IniLabels(i);
    %leafNode = find_nearest_leaf_node(gW,gW.X(i,:)); %
    j_max = gMRA.Scales(leafNode);
    
    if isfield(imgOpts, 'X0') && imgOpts.isCompressed
        nSubfigs = j_max+2; % number of subfigures to display
    else
        nSubfigs = j_max+1;
    end
    nFigsPerRow = ceil(nSubfigs/3);
    
    figure;
    for j = 1:j_max
        DataProjection = IGWTScalej( gMRA, XGWT, j );
%         if imgOpts.isCompressed
%             X_approx = imgOpts.U(:,1:D)*DataProjection(:,i)+imgOpts.cm;
%         else
            X_approx = DataProjection(:,i);
%         end
        subplot(3,double(nFigsPerRow),double(j)); imagesc(reshape(X_approx, double(imgOpts.imR), double(imgOpts.imC)))
        set(gca, 'xTick', [], 'yTick', [])
        title(num2str(j)); colormap gray
    end
    
    % original but projected
%     if imgOpts.isCompressed
%         j = j+1;
%         X_proj = imgOpts.U(:,1:D)*X(:,i)+imgOpts.cm;
%         subplot(3,double(nFigsPerRow),double(j)); imagesc(reshape(X_proj, imgOpts.imR, imgOpts.imC))
%         title 'projection'
%         set(gca, 'xTick', [], 'yTick', [])
%         colormap gray
%         
%         % original
%         if isfield(imgOpts, 'X0')
%             j = j+1;
%             X_orig = imgOpts.X0(:,i);
%             subplot(3,double(nFigsPerRow),double(j)); imagesc(reshape(X_orig, imgOpts.imR, imgOpts.imC))
%             title 'original'
%             set(gca, 'xTick', [], 'yTick', [])
%             colormap gray
%         end
%         
%     else % not compressed, then X is the original, X0 is not relevant
        j = j+1;
        X_orig = X(:,i);
        subplot(3,double(nFigsPerRow),double(j)); imagesc(reshape(X_orig, imgOpts.imR, imgOpts.imC))
        title 'original'
        set(gca, 'xTick', [], 'yTick', [])
        colormap gray
%     end
    
    %
    chain = dpath(gMRA.cp, leafNode);
    
    figure;
    for j = 1:j_max
        if ~isempty(gMRA.WavBases{chain(j)})
            wavDict_j = gMRA.WavBases{chain(j)}';
            if imgOpts.isCompressed
%                 wavDict_j = imgOpts.U(:,1:D)*wavDict_j;
                wavDict_j = gMRA.opts.Proj'*wavDict_j;
            end
            appWavDict_j = [wavDict_j; min(min(wavDict_j))*ones(imgOpts.imR, size(wavDict_j,2))];
            matWavDict_j = reshape(appWavDict_j, imgOpts.imR, []);
            matWavDict_j = matWavDict_j(:,1:end-1);
            subplot(3,ceil(length(chain)/3),double(j));
            imagesc(matWavDict_j);
            title(num2str(j)); colormap gray
            set(gca, 'xTick', [], 'yTick', [])
            %colormap(map2);
            %balanceColor;
        end
    end
    
end

return;