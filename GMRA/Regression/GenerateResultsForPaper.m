%
% Script for reproducing examples in paper
%
%% Go parallel
gcp;

addpath('./MultivariatePolynomialRegression');
addpath('./hat');

%% Set parameters
MAX_SCALES                = 100;
MSE_MC_OVERSAMPLINGFACTOR = 1;
OVERWRITE                 = true;

pExpNames           = {'SwissRoll','S-Manifold'};
pMachineNames       = {'','',''};
pPartitionMethods   = 1;

ParamsY{1}.Name     = 'sigmaY';
ParamsY{1}.Values   = num2cell([0,0.1,0.5,0.1]);
ParamsY{2}.Name     = 'Yfcn';
ParamsY{2}.Values   = {@fcn_pole,@fcn_oscpole,@fcn_randomHyperplane,@fcn_coordAxisHyperplane,@fcn_BallIndicator,@fcn_Gaussian};

ParamYList          = GenerateStructuresWithVariedParameters( ParamsY );

[~,localmachinename] = system('hostname');

pExpName = {};
for k = 1:min(length(pMachineNames),length(pExpNames));
    if isempty(pMachineNames{k}) || ~isempty(strfind(localmachinename,pMachineNames{k}))
        pExpName{end+1} = pExpNames{k};
    end
end

clear ApproxError_* Quantiles_*

for p = 1:length(pExpName)
    fprintf('\n\n --------------------\n Running data set %s \n ------------------',pExpName{p});
    
    pMainDir = '.';
    
    switch( lower(pExpName{p}) )
        case 'swissroll'
            Params{1}.Name   = 'd';
            Params{1}.Values = num2cell([2]);
            Params{2}.Name   = 'n';
            Params{2}.Values = num2cell([8000,16000,32000,64000,128000,256000,512000]);
            Params{3}.Name   = 'sigmaX';
            Params{3}.Values = num2cell([0,0.01,0.1]);                                  % will be divided by \sqrt{D} later
            
            GenerateDataSetsName = 'SwissRoll';
        case 's-manifold'
            Params{1}.Name   = 'd';
            Params{1}.Values = num2cell([2,4,6,8]);
            Params{2}.Name   = 'n';
            Params{2}.Values = num2cell([8000,16000,32000,64000,128000,256000,512000]);
            Params{3}.Name   = 'sigmaX';
            Params{3}.Values = num2cell([0,0.01,0.1]);                                  % will be divided by \sqrt{D} later
            
            GenerateDataSetsName = 'S-Manifold';
        otherwise
            return;
    end
    
    ParamList = GenerateStructuresWithVariedParameters( Params );
    
    if exist('outputs')~=9, mkdir('outputs');   end
    save([pMainDir '/outputs/' pExpName{p} '_ParamList'],'ParamList','Params');
       
    GMRAopts                        = cell(1,length(ParamList));
    Timings(1,length(ParamList))    = struct('GMRA',[],'Regression',[]);
    
    lExpName = pExpName{p};
    
    %% Run data sets with the various parameters
    for i = 1:length(ParamList),
        fprintf('\n Data Set %d/%d: %s ------------------------------------------',i,length(ParamList),pExpName{p});
        
        if (~OVERWRITE) && (exist(sprintf('%s/outputs/%s_X%.2d_idxs.mat',pMainDir,lExpName,i)))
            fprintf('\n Skipping %s_X%.2d...',lExpName,i);
            continue;
        end
        
        ParamList{i},
        %% Generate data
        fprintf('\n\t Constructing data set ...');
        [X_train, GMRAopts, imgOpts, ~,Labels]  = GenerateData_and_SetParameters(pExpName{p},ParamList{i}.n,ParamList{i}.d);        
        [X_test, ~, imgOpts, ~,Labels]          = GenerateData_and_SetParameters(pExpName{p},ParamList{i}.n,ParamList{i}.d);        
        fprintf('done.');
        
        X_train     = single(X_train);
        X_test      = single(X_test);
        
        % Add noise
        X0_train    = X_train;
        X0_test     = X_test;
        X_train     = X_train+(ParamList{i}.sigmaX/sqrt(size(X_train,1)))*randn(size(X_train));
        X0_test     = X_test+(ParamList{i}.sigmaX/sqrt(size(X_test,1)))*randn(size(X_test));

        %% Construct GMRA
        fprintf('\n\t Constructing GMRA ...');
        GMRAopts.GWTversion      = 0;                                                                                           % Set GMRA options
        if pPartitionMethods == 0
            GMRAopts.PartitionType      = 'nesdis';
            GMRAopts.smallestMetisNet   = 10;
        else
            GMRAopts.PartitionType              = 'covertree';
            GMRAopts.CoverTreeExpandOpts        = struct('ExtRange','max');
            if isfield(GMRAopts,'ManifoldDimension') && (GMRAopts.ManifoldDimension>0)
                if GMRAopts.ManifoldDimension==1,
                    theta = 0.75;
                else
                    theta = 1-1/(2*GMRAopts.ManifoldDimension);
                end
                GMRAopts.CoverTreeBuildOpts      = struct( 'theta'     , theta, ...
                    'numlevels' , max([1,int32(round(log(size(X_train,2)/(1+GMRAopts.ManifoldDimension))/log(1/theta)))]), ...
                    'minlevel'  , int32(0), 'NTHREADS'  , int32(feature('numcores')), 'BLOCKSIZE' , int32(2048));
                GMRAopts.CoverTreeTrimOpts       = struct( 'TrimType','Size','TrimValue',int32(GMRAopts.ManifoldDimension));
            else
                theta = 0.98;
                GMRAopts.CoverTreeBuildOpts      = struct( 'theta',theta,'numlevels',2*max([1,int32(round(0.5*log(1+size(X_train,2))/log(1/theta)))]),'minlevel',int32(0), ...
                    'NTHREADS',int32(feature('numcores')),'BLOCKSIZE',int32(2048));
                GMRAopts.CoverTreeTrimOpts       = struct( 'TrimType','Size','TrimValue',int32(min([size(X_train,2),10])));
            end
            GMRAopts.CoverTreeOpts.RefineCoverTree = false;
            GMRAopts.CoverTreeOpts.MaxSamples = Inf;
        end
        
        GMRAopts.parallel                   = false;
        GMRAopts.ComputeWavelets            = true;
        GMRAopts.ConstructGraph             = false;
        GMRAopts.threshold0                 = 0.8;
        GMRAopts.threshold1                 = 1e-6;
        GMRAopts.addTangentialCorrections   = true;
        GMRAopts.precision                  = 0.1;
        % Parameters for diffusion maps
        GMRAopts.ConstructGraph             = false;                                                                            % Change this to construct graph and diffusion embedding
        GMRAopts.knn                        = 50;
        GMRAopts.knnAutotune                = 20;
        GMRAopts.graphNormalization         = 'beltrami';
        GMRAopts.addTangentialCorrections   = false;
        GMRAopts.ComputeWavelets            = false;
        
        if false,
            GMRAopts.CoverTreeBuildOpts.numlevels = 100;                                                                        % Go very deep with covertree and GMRA
            GMRAopts.CoverTreeTrimOpts = 1;
        end

        Timings(i).GMRA         = cputime;
        gMRA                    = GMRA( X_train, GMRAopts );
        Timings(i).GMRA         = cputime-Timings(i).GMRA;
        fprintf('done.');   
                
        %% Generate functions Y
        Y0_train    = zeros(size(X_train,2),length(ParamYList));
        Y_train     = zeros(size(X_train,2),length(ParamYList));
        Y0_test     = zeros(size(X_train,2),length(ParamYList));
        Y_test      = zeros(size(X_train,2),length(ParamYList));
        for l = 1:length(ParamYList)
            %eval(sprintf('Yfcn=@%s;',ParamYList{l}.Yfcn));
            Yfcn            = ParamYList{l}.Yfcn;
            Y0_train(:,l)   = Yfcn( X_train );
            Y_train(:,l)    = Y0_train(:,l)+ParamYList{l}.sigmaY*randn(size(Y_train,1),1);
            Y_test(:,l)     = Yfcn( X_test );
        end
                
        %% Perform GMRA regression
        fprintf('\n\t Estimating %d regression functions...',length(ParamYList));tic        
        GMRARegressionOpts                      = struct('gMRA',gMRA,'GMRAopts',GMRAopts,'degree',[0,1,2],'noProj',false);
        Y_GMRA                                  = GMRARegression(X_train, Y_train, GMRARegressionOpts );
        fprintf('done. (%.4f sec)',toc);
        
        %% Do predictions with regression on uniform partitions
        fprintf('\n\t Evaluating regression estimator at test points...');tic
        GMRARegressionTest_opts                 = struct('allscales',true);
        Yhat                                    = GMRARegressionTest( Y_GMRA, X_test, GMRARegressionTest_opts );
        fprintf('done. (%.4f sec)',toc);
        
        %% Do predictions with regression on adaptive partitions
        fprintf('\n\t Evaluating adaptive regression estimator at test points...');tic
        GMRARegressionTest_opts                 = struct('allscales',false,'kappa',[],'noProj',GMRARegressionOpts.noProj);
        GMRARegressionTest_opts.XGWT            = Yhat.XGWT;
        GMRARegressionTest_opts.Scores          = Yhat.Scores;
        GMRARegressionTest_opts.ScalCoeffsExt   = Yhat.ScalCoeffsExt;
        GMRARegressionTest_opts.Ypreds          = Yhat.Ypreds;
        Yhat_adapt                              = GMRARegressionTest( Y_GMRA, X_test, GMRARegressionTest_opts );
        fprintf('done. (%.4f sec)',toc);
        
        %% Compute prediction error
        unifError                   = ComputeRegressionError(Y_test,Yhat);
        adaptError                  = ComputeRegressionError(Y_test,Yhat_adapt);
        Timings(i).Regression       = toc;
        Timings(i).Regression       = Timings(i).Regression/length(ParamYList);
        
        Info(i).ExpName             = pExpName{p};
        Info(i).ParamList           = ParamList{i};
        Info(i).GMRAopts            = GMRAopts;
        Info(i).GMRARegressionOpts  = GMRARegressionOpts;
        Info(i).Error.unifError     = unifError;
        Info(i).Error.adaptError    = adaptError;
        Info(i).Yhat                = Yhat;
        Info(i).Yhat_adapt          = Yhat_adapt;
        Info(i).Y0_test             = Y0_test;
        Info(i).Y_test              = Y_test;
        Info(i).Yhat_adapt.XGWT     = [];                                                                                       % This is a duplicate, same as Info(i).Yhat.XGWT
        save([pMainDir '/outputs/' pExpName{p} '_Results'],'Info','-v7.3');
    end
end
save([pMainDir '/outputs/' pExpName{p} '_Timings'],'Timings');

