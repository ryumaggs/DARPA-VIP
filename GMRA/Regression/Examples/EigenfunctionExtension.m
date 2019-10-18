%
% Script for testing GMRA regression vs. Nystrom for eigenfunction extension
%
DEG = 2;

COMPAREWITHNOPROJ           = false;
COMPARE_NYSTROMRIDGE        = false;
COMPARE_NYSTROMEIGEN        = false;
COMPARE_RSVM                = false;
COMPARE_NN                  = true;
COMPARE_RTREES              = true;

DISPLAY                     = false;
DISPLAY_ADAPTIVE_PARTITIONS = false;

N_EIGENFUNCTIONS            = 100;

%% Setup
SelectDataSetAndOptions                                                                                                         %% Select data set and other options
GenerateTrainAndTestData                                                                                                        %% Generate train and test samples

sigmaY = input('\n Noise level on training values [default=0]:          ');
if isempty(sigmaY), sigmaY = 0; end


%% Set GMRA parameters
SetGMRAParameters
if ~isfield(GMRAopts,'ManifoldDimension') || GMRAopts.ManifoldDimension==0, GMRAopts.ManifoldDimension = 0; end

if true
    GMRAopts.CoverTreeBuildOpts.numlevels = int32(100);                                                                                % Go very deep with covertree and GMRA
    GMRAopts.CoverTreeTrimOpts.Size = 5;
end

fprintf('\n The data set X consists of n=%d points in D=%d dimensions, of which %d%% is test.',size(X,2),size(X,1),Params.TestPercent);

%% Construct the Geometric Multi-Resolution Analysis
fprintf('\n Constructing GMRA on train data...');tic
GMRAopts.addTangentialCorrections   = false;
GMRAopts.ComputeWavelets            = false;
gMRA                                = GMRA( X_train, GMRAopts );
fprintf('\b done. (%.4f sec)',toc);

%% Compute Diffusion Map on all data, and on training data
GraphDiffOpts               = OutParams.GenerateDataSets.GraphDiffOpts;
GraphDiffOpts.kEigenVecs    = N_EIGENFUNCTIONS;
G                           = GraphDiffusion( X,0, GraphDiffOpts );
GraphDiffOptsTrain          = GraphDiffOpts;
GraphDiffOptsTrain.DistInfo.CoverTree = gMRA.CoverTree;
G_train                     = GraphDiffusion( X_train,0, GraphDiffOptsTrain );

%% Prepare the functions Y = eigenfunctions
Y       = G.EigenVecs;
Y_traino = Y(TrainAndTestOpts.TrainIdxs,:);
if sigmaY>0
    Y_train = Y_traino+sigmaY*randn(size(Y_traino));
else
    Y_train = Y_traino;
end
Y_test  = Y(TrainAndTestOpts.TestIdxs,:);

%% Perform GMRA regression
fprintf('\n Estimating regression function...');tic
GMRARegressionOpts              = struct('gMRA',gMRA,'GMRAopts',GMRAopts,'degree',[0:DEG],'noProj',false);
if exist('Labels','var') && ~isempty(Labels) && length(unique(Labels))<20                                                       % MM: This is a hack currently for classification...
    GMRARegressionOpts.isClassificationProblem = true;
else
    GMRARegressionOpts.isClassificationProblem = false;
end
gMRAReg                                 = GMRARegression(X_train, Y_train, GMRARegressionOpts );
fprintf('done. (%.4f sec)',toc);

%% Do predictions with regression on uniform partitions
fprintf('\n Evaluating uniform regression estimator at test points...');tic
GMRARegressionTest_opts                 = struct('allscales',true); 
Yhat                                    = GMRARegressionTest( gMRAReg, X_test, GMRARegressionTest_opts );
fprintf('done. (%.4f sec)',toc);

%% Do predictions with regression on adaptive partitions
fprintf('\n Evaluating adaptive regression estimator at test points...');tic
GMRARegressionTest_opts                 = struct('allscales',false,'kappa',[],'noProj',GMRARegressionOpts.noProj);
GMRARegressionTest_opts.XGWT            = Yhat.XGWT;
GMRARegressionTest_opts.Scores          = Yhat.Scores;
GMRARegressionTest_opts.ScalCoeffsExt   = Yhat.ScalCoeffsExt;
GMRARegressionTest_opts.Ypreds          = Yhat.Ypreds;
Yhat_adapt                              = GMRARegressionTest( gMRAReg, X_test, GMRARegressionTest_opts );
fprintf('done. (%.4f sec) \n',toc);

if COMPAREWITHNOPROJ
    fprintf('\n Constructing regression without local projections...');
    GMRARegressionOptsNoProj                        = struct('gMRA',gMRA,'GMRAopts',GMRAopts,'degree',[0,1,2],'noProj',false);
    GMRARegressionOptsNoProj.noProj                 = true;
    gMRARegNoProj                                   = GMRARegression(X_train, Y_train, GMRARegressionOptsNoProj );
    GMRARegressionTest_optsNoProj                   = struct('allscales',true);
    YhatNoProj                                      = GMRARegressionTest( gMRARegNoProj, X_test, GMRARegressionTest_optsNoProj );
    GMRARegressionTest_optsNoProj                   = struct('allscales',false,'kappa',[],'noProj',GMRARegressionOptsNoProj.noProj);
    GMRARegressionTest_optsNoProj.XGWT              = YhatNoProj.XGWT;
    GMRARegressionTest_optsNoProj.Scores            = YhatNoProj.Scores;
    GMRARegressionTest_optsNoProj.ScalCoeffsExt     = YhatNoProj.ScalCoeffsExt;
    GMRARegressionTest_optsNoProj.Ypreds            = YhatNoProj.Ypreds;
    Yhat_adaptNoProj                                = GMRARegressionTest( gMRARegNoProj, X_test, GMRARegressionTest_optsNoProj );
    fprintf('\n done. \n');
end


%% Test correctedness of adaptive partitions
if false
    % Check that there are no overlapping elements in the partition
    for l = 1:length(Yhat_adapt.partition)
        for zk = 1:length(Yhat_adapt.partition{l})
            zpathtoroot=dpath(gMRA.cp,Yhat_adapt.partition{l}(zk));
            zidxs=intersect(zpathtoroot(1:end-1),Yhat_adapt.partition{l});
            if ~isempty(zidxs), fprintf('\n\tOverlapping elements in partition %d:',l);zidxs, end
        end
    end
    % Check that the whole space is covered
    for l = 1:length(Yhat_adapt.partition)
        for zk = 1:length(gMRA.LeafNodes)
            zpathtoroot=dpath(gMRA.cp,gMRA.LeafNodes(zk));
            zidxs=intersect(zpathtoroot,Yhat_adapt.partition{l});
            if isempty(zidxs)
                fprintf('\n\tLeaf %d not covered',zk);
            elseif length(zidxs)>1
                fprintf('\n\tLeaf %d covered multiple times',zk);
            end
        end
    end
end

%% Compute prediction error
if ~gMRAReg.RegressionOpts.isClassificationProblem
    unifError   = ComputeRegressionError(Y_test,Yhat);
    adaptError  = ComputeRegressionError(Y_test,Yhat_adapt);
else                                                                                                                            % It's probably a classification problem
    unifError   = ComputeRegressionError(Y_test,Yhat,unique(Y));
    adaptError  = ComputeRegressionError(Y_test,Yhat_adapt,unique(Y));
end
if COMPAREWITHNOPROJ
    if ~gMRAReg.RegressionOpts.isClassificationProblem
        unifErrorNoProj   = ComputeRegressionError(Y_test,YhatNoProj);
        adaptErrorNoProj  = ComputeRegressionError(Y_test,Yhat_adaptNoProj);
    else                                                                                                                            % It's probably a classification problem
        unifErrorNoProj   = ComputeRegressionError(Y_test,YhatNoProj,unique(Y));
        adaptErrorNoProj  = ComputeRegressionError(Y_test,Yhat_adaptNoProj,unique(Y));
    end
end


%% Visualizes the best partitions
if DISPLAY_ADAPTIVE_PARTITIONS
    [~,minidx] = min(adaptError.rel(1,:));
    GMRA_VisualizePartition(X_train,gMRA,Yhat_adapt.partition{2}(1:10));
end

drawnow;
%[zzX,zzY,zzT,zzAUC] = perfcurve(Y_test(:,1),Yhat.Yhat{1}{end}(:,1)/max(Yhat.Yhat{1}{end}(:,1)),1);figure;plot(zzX,zzY);

%% Compare with regression trees
clear Y_test_regtree Error_regTree RegTree
if COMPARE_RTREES
    fprintf('\n Constructing regression trees...');
    for k = size(Y,2):-1:1
        if ~gMRAReg.RegressionOpts.isClassificationProblem
            RegTree{k} = fitrtree(X_train',Y_train(:,k),'CrossVal','off');
            Y_test_regtree(:,k) = predict(RegTree{k},X_test');
            Error_regTree(k) = norm(Y_test(:,k)-Y_test_regtree(:,k))/norm(Y_test(:,k));
        else
            RegTree{k} = fitctree(X_train',Y_train(:,k),'CrossVal','off');
            Y_test_regtree(:,k) = predict(RegTree{k},X_test');
            Error_regTree(k) = sum(Y_test(:,k)~=Y_test_regtree(:,k))/size(Y_test,1);
        end
    end
    fprintf('done.\n');
end

%% Compare with NN
clear Y_testNN Error_NN
if COMPARE_NN
    fprintf('\n Performing NN regression...');
    for k = size(Y,2):-1:1
        if isfield(gMRA,'CoverTree')
            Y_testNN{k} = NNclassify( X_train, Y_train(:,k), X_test, 1:20, gMRA.CoverTree, gMRAReg.RegressionOpts.isClassificationProblem );
        else
            Y_testNN{k} = NNclassify( X_train, Y_train(:,k), X_test, 1:20, [], gMRAReg.RegressionOpts.isClassificationProblem );
        end
        if gMRAReg.RegressionOpts.isClassificationProblem
            for l = size(Y_testNN{k},2):-1:1
                Error_NN(l,k) = sum(Y_testNN{k}(:,l)~=Y_test(:,k))/length(Y_test(:,k));
            end
        else
            for l = size(Y_testNN{k},2):-1:1
                Error_NN(l,k) = norm(Y_testNN{k}(:,l)-Y_test(:,k))/norm(Y_test(:,k));
            end
        end
    end
    fprintf('done.\n');
end


%% Compare with Ridge regression a la Nystrom
if COMPARE_NYSTROMRIDGE
    fprintf('\n Performing ridge regression...');
    clear nystromCoRe_config nystromCoRe_training_output nystromCoRe_prediction_output Error_nystromCoRe
    [medianmindist,medianmediandist,minmindist] = EstimateMinPairwiseDistance( X );
    nystrcomCoRe_kernelParameter                = linspace(minmindist,medianmediandist/2,5);
    for l = 1:length(nystrcomCoRe_kernelParameter)
        nystromCoRe_config(l)                   = config_set( 'recompute', 1, 'kernel.kernelParameter',nystrcomCoRe_kernelParameter(l) );
        for k = size(Y,2):-1:1
            nystromCoRe_training_output(k,l)    = nystromCoRe_train ( X_train' , Y_train(:,k), nystromCoRe_config(l) );
            nystromCoRe_prediction_output(k,l)  = nystromCoRe_test  ( X_test' , Y_test(:,k) , nystromCoRe_training_output(k,l) );
            Error_nystromCoRe(l,k)              = norm(nystromCoRe_prediction_output(k,l).YtePred-Y_test(:,k))/norm(Y_test(:,k));
        end
    end
    fprintf('done.');
    %figure;plot(nystrcomCoRe_kernelParameter,log10(Error_nystromCoRe'));axis tight;xlabel('\sigma');
end

%% Compare with Nystrom method for extension of eigenfunctions of the Laplacian
if COMPARE_NYSTROMEIGEN
    Y_NystromEigen      = Nystrom(X_train,X_test,G_train,GraphDiffOpts.kNN);
    Error_NystromEigen  = min([sqrt(sum((Y_test-Y_NystromEigen).^2,1))./sqrt(sum(Y_test.^2,1));sqrt(sum((Y_test+Y_NystromEigen).^2,1))./sqrt(sum(Y_test.^2,1))]);
end

%% Now compare with regression SVM
if COMPARE_RSVM
    fprintf('\n Constructing regression SVM''s...');
    clear ParamsRSVM ParamRSVMList RSVM Y_test_RSVM Error_RSVM
    ParamsRSVM{1}.Name     = 'KernelFunction';
    ParamsRSVM{1}.Values   = {'linear','gaussian','rbf','polynomial'};
    ParamsRSVM{2}.Name     = 'PolynomialOrder';
    ParamsRSVM{2}.Values   = {1};
    ParamRSVMList          = GenerateStructuresWithVariedParameters( ParamsRSVM );
    L = length(ParamRSVMList);
    for l = 1:L
        if strcmp(ParamRSVMList{l}.KernelFunction,'polynomial')
            for punif = 2:3
                ParamRSVMList{end+1} = ParamRSVMList{l};
                ParamRSVMList{end}.PolynomialOrder = punif;
            end
        end
    end
    
    for l = 1:length(ParamRSVMList)
        for k = size(Y,2):-1:1
            if strcmp(ParamRSVMList{l}.KernelFunction,'polynomial')
                RSVM{l,k} = fitrsvm(X_train',Y_train(:,k),'KernelFunction',ParamRSVMList{l}.KernelFunction,'BoxConstraint',max(abs(Y_train(:,k))),'PolynomialOrder',ParamRSVMList{l}.PolynomialOrder,'CrossVal','off');
            else
                RSVM{l,k} = fitrsvm(X_train',Y_train(:,k),'KernelFunction',ParamRSVMList{l}.KernelFunction,'BoxConstraint',max(abs(Y_train(:,k))),'CrossVal','off');
            end
            Y_test_RSVM(:,l,k) = predict(RSVM{l,k},X_test');
            Error_RSVM(l,k) = norm(Y_test(:,k)-squeeze(Y_test_RSVM(:,l,k)))/norm(Y_test(:,k));
        end
    end
    fprintf('done.');
    Error_RSVM
end

%% Numerical results
for k = 1:size(Y_test,2)
    fprintf('\n\n SUMMARY of performance for function %d --------------------------------------------------------',k)
    fprintf('\n Best result from multiscale uniform estimators of degree %d: %f',[0:DEG;min(squeeze(unifError.rel(:,:,k)),[],2)'])
    fprintf('\n Best result from multiscale adaptive estimators of degree %d: %f',[0:DEG;min(squeeze(adaptError.rel(:,:,k)),[],2)'])
    if COMPAREWITHNOPROJ
        fprintf('\n Best result from multiscale uniform non-projected estimators of degree %d: %f',[0:DEG;min(squeeze(unifErrorNoProj.rel(:,:,k)),[],2)'])
        fprintf('\n Best result from multiscale adaptive non-projected estimators of degree %d: %f',[0:DEG;min(squeeze(adaptErrorNoProj.rel(:,:,k)),[],2)'])
    end
    if COMPARE_RTREES
        fprintf('\n Best result from regression and classification trees: %f',Error_regTree(k));
    end
    if COMPARE_NN
        fprintf('\n Best result from NN estimators: %f',min(squeeze(Error_NN(:,k))));
    end
    if COMPARE_NYSTROMRIDGE
        fprintf('\n Best result from NystromCoRe: %f',min(squeeze(Error_nystromCoRe(:,k))));
    end
    if COMPARE_NYSTROMEIGEN
        fprintf('\n Best result from NystromEigen: %f',Error_NystromEigen(k));
    end
end

%% Plot results
clear min_by_*
for deg = DEG:-1:0
    [min_by_deg_unif(deg+1,:),min_by_deg_idx_unif(deg+1,:)]     = min(squeeze(unifError.rel(deg+1,:,:)),[],1);
    min_by_deg_partitionsize_unif(deg+1,:)                      = cellfun(@length,Yhat.partition{deg+1}(min_by_deg_idx_unif(deg+1,:)));    
    [min_by_deg_adapt(deg+1,:),min_by_deg_idx_adapt(deg+1,:)]   = min(squeeze(adaptError.rel(deg+1,:,:)),[],1);
    for i = 1:size(Y_test,2)
        min_by_deg_partitionsize_adapt(deg+1,i)                 = length(Yhat_adapt.partition{deg+1}{min_by_deg_idx_adapt(deg+1,i),i});
        min_Y_unif{deg+1}(:,i)                                  = Yhat.Yhat{deg+1}{min_by_deg_idx_unif(deg+1,i)}(:,i);
        min_Y_adapt{deg+1}(:,i)                                 = Yhat_adapt.Yhat{deg+1}{min_by_deg_idx_unif(deg+1,i)}(:,i);
    end
end
% Plot optimal error
bigFig;
punif=plot(1:N_EIGENFUNCTIONS,min_by_deg_unif);set(punif(1),'Color',[0,0,1],'Marker','.');set(punif(2),'Color',[0,0,0.8],'Marker','o');set(punif(3),'Color',[0,0,0.6],'Marker','x');hold on;
padapt=plot(1:N_EIGENFUNCTIONS,min_by_deg_adapt);set(padapt(1),'Color',[1,0,0],'Marker','.');set(padapt(2),'Color',[0.8,0,0],'Marker','o');set(padapt(3),'Color',[0.6,0,0],'Marker','x');
plot_list = [punif;padapt];
legendstrings = {'unif_0','unif_1','unif_2','adapt_0','adapt_1','adapt_2'};
if COMPARE_RTREES
    prtrees = plot(1:N_EIGENFUNCTIONS,Error_regTree);set(prtrees,'Color',[0,1,0],'Marker','.');
    plot_list = [plot_list;prtrees];
    legendstrings{end+1} = 'R-Trees';
end
if COMPARE_NN
    pnns = plot(1:N_EIGENFUNCTIONS,min(Error_NN,[],1));set(pnns,'Color',[0,1,0],'Marker','o');
    plot_list = [plot_list;pnns];
    legendstrings{end+1} = 'k-NN';
end
if COMPARE_RSVM
    prsvm = plot(1:N_EIGENFUNCTIONS,min(Error_RSVM,[],1));set(prsvm,'Color',[0,1,0],'Marker','x');
    plot_list = [plot_list;prsvm];
    legendstrings{end+1} = 'R-SVM';
end
if COMPARE_NYSTROMRIDGE
    pnystromridge = plot(1:N_EIGENFUNCTIONS,min(Error_RSVM,[],1));set(pnystromridge,'Color',[0,1,0],'Marker','*');
    plot_list = [plot_list;pnystromridge];
    legendstrings{end+1} = 'Nystrom Ridge';
end

peigs = plot(1:N_EIGENFUNCTIONS,1-G.EigenVals);
set(peigs,'Color',[0,0,0],'LineStyle','--');
plot_list = [plot_list;peigs];
legendstrings{end+1} = '\lambda';

legend(plot_list,legendstrings,'Location','NorthWest');
xlabel('Eigenfunction index'); ylabel('Average L^2 normalized error');
title('Performance of various algorithms as a function of eigenfunction index');
axis tight
% Size of optimal (over choice of scale/kappa) partitions
%figure;subplot(1,2,1);imagesc(min_by_deg_partitionsize_unif);colorbar;subplot(1,2,2);imagesc(min_by_deg_partitionsize_adapt);colorbar'

%% Matrix of inner products between true and estimate eigenvectors
bigFig;
for deg = 1:DEG+1
    subplot(2,DEG+1,deg); imagesc(log10(abs(Y_test'*min_Y_unif{deg})));colorbar; title(sprintf('Unif_{%d}',deg-1));xlabel('Y_{pred}');ylabel('Y');
    subplot(2,DEG+1,deg+DEG+1); imagesc(log10(abs(Y_test'*min_Y_adapt{deg})));colorbar; title(sprintf('Adapt_{%d}',deg-1));xlabel('Y_{pred}');ylabel('Y');
end

fprintf('\n');