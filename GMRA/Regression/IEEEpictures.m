%
% Reduced/modified version of RunExamples
% to generate pictures for IEEE paper
%

COMPAREWITHNOPROJ = true;

%% Setup
SelectDataSetAndOptions                                                                                                         %% Select data set and other options
GenerateTrainAndTestData                                                                                                        %% Generate train and test samples
%% Set GMRA parameters
SetGMRAParameters
if ~isfield(GMRAopts,'ManifoldDimension') || GMRAopts.ManifoldDimension==0, GMRAopts.ManifoldDimension = 10; end

if true,
    GMRAopts.CoverTreeBuildOpts.numlevels = 100;                                                                                % Go very deep with covertree and GMRA
    GMRAopts.CoverTreeTrimOpts.Size = 5;
end


fprintf('\n The data set X consists of n=%d points in D=%d dimensions, of which %d%% is test.',size(X,2),size(X,1),pTestPercent);

%% Construct the Geometric Multi-Resolution Analysis
fprintf('\n Constructing GMRA on train data...');tic
GMRAopts.addTangentialCorrections   = false;
GMRAopts.ComputeWavelets            = false;
%X_train = X_train + 0.01*randn(size(X_train));
%X_test = X_test + 0.01*randn(size(X_test));
gMRA    = GMRA( X_train, GMRAopts );
fprintf('\b done. (%.4f sec)',toc);

%% Generate the function Y
if ~exist('Labels','var') || isempty(Labels)
    %    sigmaY = input('\n Noise level on Y [default=0]:          ');
    %    if isempty(sigmaY), sigmaY = 0; end
    sigmaY = 0.05;
    Yfcn    = {@fcn_norm_withspike};%,@fcn_BallIndicator,@fcn_oscpole,@fcn_norm_withbigspike};
    Y       = [];
    Yo      = [];
    for l = 1:length(Yfcn)
        Ynew    = Yfcn{l}( X );
        Yo      = [Yo,Ynew];
        Y       = [Y,Ynew+sigmaY*randn(size(Ynew))];
    end
    clear Ynew;
else
    sigmaY = 0;
    [uniqueLabels,~,Y] = unique(Labels);
    Yo = Y;
end
Y_train = Y(TrainAndTestOpts.TrainIdxs,:);
Y_test  = Yo(TrainAndTestOpts.TestIdxs,:);

%% Perform GMRA regression
fprintf('\n Estimating regression function...');tic
GMRARegressionOpts              = struct('gMRA',gMRA,'GMRAopts',GMRAopts,'degree',[0,1,2],'noProj',false);
if exist('Labels','var') && ~isempty(Labels) && length(unique(Labels))<20                                                       % MM: This is a hack currently for classification...
    GMRARegressionOpts.isClassificationProblem = true;
else
    GMRARegressionOpts.isClassificationProblem = false;
end
gMRAReg                                 = GMRARegression(X_train, Y_train, GMRARegressionOpts );
fprintf('done. (%.4f sec)',toc);

%% Do predictions with regression on uniform partitions
fprintf('\n Evaluating regression estimator at test points...');tic
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

fprintf('done. (%.4f sec)',toc);

%% Test correctedness of adaptive partitions
if false,
    % Check that there are no overlapping elements in the partition
    for l = 1:length(Yhat_adapt.partition)
        for zk = 1:length(Yhat_adapt.partition{l});
            zpathtoroot=dpath(Y_GMRA.gMRA.cp,Yhat_adapt.partition{l}(zk));
            zidxs=intersect(zpathtoroot(1:end-1),Yhat_adapt.partition{l});
            if ~isempty(zidxs), fprintf('\n\tOverlapping elements in partition %d:',l);zidxs, end
        end
    end
    % Check that the whole space is covered
    for l = 1:length(Yhat_adapt.partition)
        for zk = 1:length(Y_GMRA.gMRA.LeafNodes)
            zpathtoroot=dpath(Y_GMRA.gMRA.cp,Y_GMRA.gMRA.LeafNodes(zk));
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
% %% Visualizes the best partitions
% if false
%     [~,minidx] = min(adaptError.rel(1,:));
%     GMRA_VisualizePartition(X_train,gMRA,Yhat_adapt.partition{2}(1:10));
% end

drawnow;

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
    if ~gMRAReg.RegressionOpts.isClassificationProblem
        unifErrorNoProj   = ComputeRegressionError(Y_test,YhatNoProj);
        adaptErrorNoProj  = ComputeRegressionError(Y_test,Yhat_adaptNoProj);
    else                                                                                                                            % It's probably a classification problem
        unifErrorNoProj   = ComputeRegressionError(Y_test,YhatNoProj,unique(Y));
        adaptErrorNoProj  = ComputeRegressionError(Y_test,Yhat_adaptNoProj,unique(Y));
    end    
    fprintf('\n done.');
end



%% Compare with regression trees
fprintf('\n Constructing regression trees...');
clear Y_test_regtree Error_regTree RegTree
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

%% Compare with NN
clear Y_testNN Error_NN
for k = size(Y,2):-1:1,
    Y_testNN{k} = NNclassify( X_train, Y_train(:,k), X_test, 1:20, gMRA.CoverTree, gMRAReg.RegressionOpts.isClassificationProblem );
    if gMRAReg.RegressionOpts.isClassificationProblem
        for l = size(Y_testNN{k},2):-1,1
            Error_NN(k,l) = sum(Y_testNN{k}(:,l)~=Y_test(:,k))/length(Y_test(:,k));
        end
    else
        for l = size(Y_testNN,2):-1:1
            Error_NN(k,l) = norm(Y_testNN{k}(:,l)-Y_test(:,k))/norm(Y_test(:,k));
        end
    end
end


%% Compare with Ridge regression a la Nystrom
fprintf('\n Performing ridge regression...');
clear nystromCoRe_config nystromCoRe_training_output nystromCoRe_prediction_output Error_nystromCoRe
%nystrcomCoRe_kernelParameter            = linspace(0.01,0.4,20);
[medianmindist,medianmediandist,minmindist] = EstimateMinPairwiseDistance( X );
nystrcomCoRe_kernelParameter                = linspace(minmindist,medianmediandist,10);
for l = 1:length(nystrcomCoRe_kernelParameter)
    nystromCoRe_config(l)                   = config_set( 'recompute', 1, 'kernel.kernelParameter',nystrcomCoRe_kernelParameter(l) );
    for k = size(Y,2):-1:1
        nystromCoRe_training_output(k,l)    = nystromCoRe_train ( X_train' , Y_train(:,k), nystromCoRe_config(l) );
        nystromCoRe_prediction_output(k,l)  = nystromCoRe_test  ( X_test' , Y_test(:,k) , nystromCoRe_training_output(k,l) );
        Error_nystromCoRe(k,l)              = norm(nystromCoRe_prediction_output(k,l).YtePred-Y_test(:,k))/norm(Y_test(:,k));
    end
end
fprintf('done.');
%figure;plot(nystrcomCoRe_kernelParameter,log10(Error_nystromCoRe));axis tight;xlabel('\sigma');


%% Now compare with regression SVM
if false
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
            for p = 2:3
                ParamRSVMList{end+1} = ParamRSVMList{l};
                ParamRSVMList{end}.PolynomialOrder = p;
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

% %% Numerical results
% for k = 1:size(Y_test,2)
%     fprintf('\n SUMMARY of performance for function %d --------------------------------------------------------',k)
%     fprintf('\n Best result from multiscale uniform estimators, of varying degrees: \n');min(squeeze(unifError.rel(:,:,k)),[],2)'
%     fprintf('\n Best result from multiscale adaptive estimators, of varying degrees: \n');min(squeeze(adaptError.rel(:,:,k)),[],2)'
%     fprintf('\n Best result from regression and classification trees: \n');Error_regTree(k),
%     fprintf('\n Best result from NN estimators: \n');min(squeeze(Error_NN(k,:))),
%     fprintf('\n Best result from NystromCoRe: \n');min(squeeze(Error_nystromCoRe(k,:))),
% end
%%
competitors = struct('alg',[],'err',[],'col',[]);
competitors.alg = {'CART','Nearest Neighbors','Nystrom','NoProj linear unif.','NoProj linear adapt.','NoProj quad. unif.','NoProj quad. adapt.'};
competitors.err = {Error_regTree,min(squeeze(Error_NN(k,:))),min(squeeze(Error_nystromCoRe(k,:))),min(unifErrorNoProj.abs(2,:)),min(adaptErrorNoProj.abs(2,:)),min(unifErrorNoProj.abs(3,:)),min(adaptErrorNoProj.abs(3,:))};
competitors.col = {'m', 'c' , 'g',[0.7,0.4,0],[0.8,0.5,0],[0.9,0.6,0],[1,0.7,0]};

%% Display accuracy
IEEEdisplay( X_test,gMRAReg, Y_test, Yhat, Yhat_adapt, unifError, adaptError, sigmaY , competitors);