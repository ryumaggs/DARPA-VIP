clear Params

%% Setup
SelectDataSetAndOptions                                                                                                         %% Select data set and other options
GenerateTrainAndTestData                                                                                                        %% Generate train and test samples
SetGMRAParameters                                                                                                               %% Set GMRA parameters

%GMRAopts.CoverTreeBuildOpts.numlevels = 100;                              % Go deep: useful for seing bias/variance trade-off
GMRAopts.CoverTreeTrimOpts.TrimValue = max(5,Params.d); 
%GMRAopts.Predictor.PredictorType = 'orthogonal';

%% GMRA construction
fprintf('\n The data set X consists of n=%d points in D=%d dimensions.',size(X,2),size(X,1));
fprintf('\n Constructing GMRA...');
gMRA = GMRA(X_train, GMRAopts);                                                                                                       %% Construct the Geometric Multi-Resolution Analysis
fprintf('done.');

%% Compute all wavelet coefficients via Fast Geometric Wavelet Transform                                                        
tic;
fprintf('\n Computing FGWT''s...');
XGWT        = FGWT(gMRA, X_train);                                                                                              % Compute FGWT of train data set
if exist('X0_test','var')
    XGWT_test   = FGWT(gMRA, X0_test);                                                                                          % Compute FGWT of test data set
end
fprintf('done. (%f secs)',toc);

%% Compute MSE for uniform and adaptive GMRA
fprintf('\n Compute MSE versus scales/partition size on traning data for uniform and adaptive GMRA...');
ErrOpts                 = struct('norm',[2,inf],'relative',true,'quantiles',[0.25,0.5,0.75]);                                   % Default option for approximation error calculation
RefinementOpts          = struct('refinementnorm','2','ifscaled',1);                
DataError_train         = AGWT_ComputeApproxErr(gMRA,XGWT,X_train,XGWT,RefinementOpts);

%% Display Results                                                                                                              % Display some results, approximation errors, etc...
fprintf('\n Compute MSE versus scales/partition size on test data for uniform and adaptive GMRA...');
%AGWT_DisplayResults( gMRA, XGWT, DataError_train );
if exist('X0_test','var')
    DataError_test      = AGWT_ComputeApproxErr(gMRA,XGWT,X0_test,XGWT_test,RefinementOpts);
end

%% Display Results                                                                                                              % Display some results, approximation errors, etc...
fprintf('\n Displaying results...');
displayopts             = struct;                                                                                               % display options                                                
displayopts.errortype   = 1;                                                                                                    % 0: absolute error 1: relative error
displayopts.errornorm   = RefinementOpts.refinementnorm;                                                                        % '2': L2 norm; 'inf': infinity norm
% AGWT_DisplayResults( gMRA, XGWT, DataError_train ,displayopts);
if exist('X0_test','var')
    AGWT_DisplayResults( gMRA, XGWT, DataError_test,displayopts);
end
if ~isfield(GMRAopts,'ManifoldDimension') || GMRAopts.ManifoldDimension==0
    dimensions=cellfun( @(x) size(x,1),gMRA.ScalBasis );
    figure;histogram(dimensions);axis tight;title('Dimensions of the Scaling Subspaces');
end
%AGWT_DisplayLinearSets(gMRA,X,XGWT_test,DataError_test);                                                                       % Display point clouds representing approximations at each scale

fprintf('\ndone.\n');

%% save data
%save([pExampleNames{pExampleIdx} 'Dim' num2str(GMRAopts.ManifoldDimension) 'S' num2str(size(X,2)) '.mat'],'-v7.3')

%% Display image patches if that's the type of data
% if isfield(gMRA.opts.GenerateDataSets,'LoadImageBatchParams')
%     figure;                                                                                                                     %% display some geometric scaling functions, at varying scales
%     for lp = 1:2
%         for j = 1:max(gMRA.Scales)
%             idxs = find(gMRA.Scales==j);
%             idxs = idxs(randperm(length(idxs),min([20,length(idxs)])));
%             maxDim = max(cellfun(@(x) size(x,1),gMRA.ScalBasis(idxs)));
%             for k = 1:length(idxs)
%                 i = randperm(size(gMRA.ScalBasis{idxs(k)},1),1);
%                 subplot(max(gMRA.Scales),20,(j-1)*20+k);
%                 imagesc(reshape(gMRA.ScalBasis{idxs(k)}(i,:),gMRA.opts.GenerateDataSets.LoadImageBatchParams.PatchSz));colormap(gray);x`
%             end
%         end
%         pause
%     end
% end


%% Compute Inverse Geometric Wavelet Transform                                                                                  % Compute IGWT (this is not used in what follows)
if false,
    fprintf('\n Computing IGWT...');
    tic;
    [Projections, tangentialCorrections] = IGWT(gMRA, XGWT);
    fprintf('done. (%f secs)',toc);
end


%% Save file fro web UI
% image parameters metadata for visualization
try
    if ~isempty(fieldnames(imgOpts))
        imgOpts.imageData               = true;
        imgOpts.Labels                  = int32(reshape(imgOpts.Labels,length(imgOpts.Labels),1));
        imgOpts.LabelNames{1}           = 'img_cat';
        imgOpts.LabelDescriptions{1}    = 'Image Label';
    else
        imgOpts.imageData = false;
    end
    
    imgOpts.LabelDataTypes{1}           = 'int32';
    imgOpts.LabelVariableTypes{1}       = 'categorical';
    imgOpts.LabelNames{1}               = 'no label';
    imgOpts.LabelDescriptions{1}        = 'no label';
    
    if ~isfield(imgOpts,'Labels') || isempty(imgOpts.Labels)
        imgOpts.Labels = ones(size(X,2),1);
    end
catch
end


%% Display diffusion embedding
if isfield(gMRA,'Graph') && ~isempty(gMRA.Graph.EigenVecs)
    figure;
    subplot(2,3,1);p=scatter3(X0(1,:),X0(2,:),X0(3,:),10,1:size(X0,2),'filled');
    subplot(2,3,2);p=scatter3(gMRA.Graph.EigenVecs(:,2),gMRA.Graph.EigenVecs(:,3),gMRA.Graph.EigenVecs(:,4),10,1:size(X,2),'filled');
    subplot(2,3,3);p=scatter3(gMRA.Graph.EigenVecs(:,3),gMRA.Graph.EigenVecs(:,4),gMRA.Graph.EigenVecs(:,5),10,1:size(X,2),'filled');
    subplot(2,3,4);p=scatter3(gMRA.Graph.EigenVecs(:,4),gMRA.Graph.EigenVecs(:,5),gMRA.Graph.EigenVecs(:,6),10,1:size(X,2),'filled');
    subplot(2,3,5);p=scatter3(gMRA.Graph.EigenVecs(:,5),gMRA.Graph.EigenVecs(:,6),gMRA.Graph.EigenVecs(:,7),10,1:size(X,2),'filled');
    subplot(2,3,6);p=scatter3(gMRA.Graph.EigenVecs(:,6),gMRA.Graph.EigenVecs(:,7),gMRA.Graph.EigenVecs(:,8),10,1:size(X,2),'filled');
end



%matlab_to_hdf5_write(gMRA, imgOpts,[pVisDir pExampleNames{pExampleIdx} '.hdf5'] );
