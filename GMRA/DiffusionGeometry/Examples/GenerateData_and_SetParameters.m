function [X, GMRAopts, imgOpts,X0,Labels,OutParams] = GenerateData_and_SetParameters(pExampleName,nSamples,d,GMRAopts,D)

X         = [];
X0        = [];
Labels    = [];
OutParams = struct();

% set GWT parameters
if nargin<4, GMRAopts = struct();   end
if nargin<5, D = max([d+10,20]);    end

% The following thresholds are used in the code construct_GMRA.m
GMRAopts.threshold1 = 1e-6;                                                                 % Threshold of singular values for determining the rank of each ( I - \Phi_{j,k} * \Phi_{j,k} ) * Phi_{j+1,k'}
GMRAopts.threshold2 = 1e-6;                                                                 % Threshold for determining the rank of intersection of ( I - \Phi_{j,k} * \Phi_{j,k} ) * Phi_{j+1,k'}

% whether to use best approximations
GMRAopts.addTangentialCorrections = true;

% whether to sparsify the scaling functions and wavelet bases
GMRAopts.sparsifying = false;
GMRAopts.sparsifying_method = 'ksvd'; % or 'spams'

% whether to split the wavelet bases into a common intersection and
% children-specific parts
GMRAopts.splitting = false;

% METIS parameters
GMRAopts.knn = 30;
GMRAopts.knnAutotune = 20;
GMRAopts.smallestMetisNet = 10;

% whether to output time
GMRAopts.verbose = 0;

% method for shrinking the wavelet coefficients
GMRAopts.shrinkage = 'hard';

% whether to merge the common part of the wavelet subspaces
% associated to the children into the scaling function of the parent.
GMRAopts.mergePsiCapIntoPhi  = false;

%% create data, and set additional parameters
imgOpts = struct();
Qrnd    = [];

if GMRAopts.verbose>0     fprintf('\nGenerating/loading %s data...', pExampleName); end
tic
switch pExampleName
    
    case 'MNIST_Digits'
        
        % generate the dataset
        dataset = struct();
        dataset.N = round(nSamples/10);
        dataset.digits = 0:9;
        dataset.projectionDimension = 1000;
        
        [X0,Labels,OutParams.GenerateDataSets] = GenerateDataSets( 'BMark_MNIST', ...
            struct('NumberOfPoints',dataset.N,'AutotuneScales',false,'MnistOpts', ...
            struct('Sampling', 'FirstN', 'QueryDigits',dataset.digits, 'ReturnForm', 'vector'))); %#ok<ASGLU>
        imgOpts.imageData       = true;
        imgOpts.imR             = OutParams.GenerateDataSets.GraphDiffOpts.imgSize(1);
        imgOpts.imC             = OutParams.GenerateDataSets.GraphDiffOpts.imgSize(2);
        imgOpts.Labels          = Labels;
        imgOpts.Coords          = OutParams.GenerateDataSets.Coords;
        
        if dataset.projectionDimension>0 && dataset.projectionDimension<imgOpts.imR*imgOpts.imC,
            if nargin<4 || ~isfield(GMRAopts,'Mean') || ~isfield(GMRAopts,'Proj') || isempty(GMRAopts.Mean) || isempty(GMRAopts.Proj)
                GMRAopts.Mean   = mean(X0,2);
                X               = bsxfun(@minus,X0,GMRAopts.Mean);
                [U,S,V]         = randPCA(X, dataset.projectionDimension);
                GMRAopts.Proj   = U';
            end
            X                   = GMRAopts.Proj*bsxfun(@minus,X0,GMRAopts.Mean);
            imgOpts.U           = GMRAopts.Proj;
            imgOpts.cm          = GMRAopts.Mean;
            imgOpts.isCompressed= true;
        else
            X = X0;
            GMRAopts.Proj = [];
            GMRAopts.Mean = [];
            imgOpts.isCompressed = false;
        end
        
        % GWT parameters that need to be set separately
        %GWTopts.ManifoldDimension = 4; % if 0, then determine locally adaptive dimensions using the following fields:
        GMRAopts.threshold0 = 0.5; % threshold for choosing pca dimension at each nonleaf node
        GMRAopts.errorType = 'relative';
        GMRAopts.precision  = .050; % only for leaf nodes
        
    case 'YaleB_Faces'
        
        load YaleB_PCA
        X = S*V'; %#ok<NODEF>
        % image parameters
        imgOpts.imageData = true;
        imgOpts.imR = 480;
        imgOpts.imC = 640;
        
        imgOpts.Labels = Labels; %#ok<NODEF>
        imgOpts.cm =  Imean;
        imgOpts.U = U; %#ok<NODEF>
        imgOpts.isCompressed = true;
        
        % GWT parameters that need to be set separately
        % GWTopts.ManifoldDimension = 4;
        GMRAopts.errorType = 'relative';
        GMRAopts.threshold0 = 0.5; % threshold for choosing pca dimension at each nonleaf node
        GMRAopts.precision  = 0.05; % only for leaf nodes
        
    case 'croppedYaleB_Faces'
        
        load extendedYaleB_crop_SVD
        dataset.projectionDimension = 500;
        %X = V(:,1:dataset.projectionDimension); %#ok<NODEF>
        X = S(1:dataset.projectionDimension,1:dataset.projectionDimension)*V(:,1:dataset.projectionDimension)'; %#ok<NODEF>
        % image parameters
        imgOpts.imageData = true;
        imgOpts.imR = 192;
        imgOpts.imC = 168;
        
        %imgOpts.Labels = Labels; %#ok<NODEF>
        imgOpts.cm =  center;
        %imgOpts.V = bsxfun(@times, U(:,1:dataset.projectionDimension), (diag(S(1:dataset.projectionDimension,1:dataset.projectionDimension)))'); %#ok<NODEF>
        imgOpts.U = U(:,1:dataset.projectionDimension); %#ok<NODEF>
        imgOpts.isCompressed = true;
        
        % GWT parameters that need to be set separately
        %GWTopts.ManifoldDimension = 4;
        GMRAopts.threshold0 = .5; % threshold for choosing pca dimension at each nonleaf node
        GMRAopts.errorType = 'relative';
        GMRAopts.precision  = .05; % only for leaf nodes
        
    case 'ScienceNews'
        
        load X20
        
        X = bsxfun(@rdivide,X,sqrt(sum(X.^2,2))); %#ok<NODEF>
        X = X';
        
        GMRAopts.ManifoldDimension = 0;
        GMRAopts.threshold0 = 0.5;
        GMRAopts.errorType = 'absolute';
        GMRAopts.precision  = 1e-5; % only for leaf nodes
        
    case 'NaturalImagePatches'
        
        load NaturalImagePatches.mat
        %X = Y(:, randsample(size(Y,2),10000)); %#ok<NODEF>
        X = single(Y);
        
        % image parameters
        imgOpts.imageData = true;
        imgOpts.imR = 16;
        imgOpts.imC = 16;
        imgOpts.isCompressed = false;
        
        % GWT parameters that need to be set separately
        %GWTopts.ManifoldDimension = 4;
        GMRAopts.errorType = 'relative';
        GMRAopts.threshold0 = 0.5; % threshold for choosing pca dimension at each nonleaf node
        GMRAopts.precision  = 0.05; % only for leaf nodes
        
    case 'IntegralOperator'
        
        [X0,Labels,OutParams.GenerateDataSets]=GenerateDataSets('Planes',   ...
            struct('NumberOfPoints',5000,'EmbedDim',3,'PlanesOpts',struct('Number',2,'Dim',[1,2], ...
            'Origins',zeros(2,3),'NumberOfPoints',[1000,4000],'Dilations',[1,5])));
        X=GenerateIntegralOperatorFromSets(X0(:,Labels==1)', X0(:,Labels==2)');
        
        %GWTopts.ManifoldDimension = 0; % if 0, then determine locally adaptive dimensions using the following fields:
        
        GMRAopts.errorType = 'relative';
        GMRAopts.threshold0 = 0.5; % threshold for choosing pca dimension at each nonleaf node
        GMRAopts.precision  = 0.01; % only for leaf nodes
        
    case 'MeyerStaircase-d'
        dataset = struct();
        dataset.name = pExampleName;
        dataset.N = nSamples;
        if nargin>=3
            dataset.d = d;
        else
            dataset.d = 2;
        end
        dataset.D = 2000;
        dataset.MeyerStepWidth=0.2;
        dataset.noiseLevel = 0/sqrt(dataset.D);
        
        % Generate data
        [X_clean,Labels,OutParams.GenerateDataSets] = GenerateDataSets( dataset.name, ...
            struct('NumberOfPoints',dataset.N,'Dim',dataset.d,'width',dataset.MeyerStepWidth, ...
            'EmbedDim',dataset.D,'NoiseType','Gaussian','NoiseParam',0) );
        
        % Add noise
        if dataset.noiseLevel>0,
            X = X_clean + dataset.noiseLevel*random('norm', 0,1, size(X_clean));
            GMRAopts.X_clean = X_clean;
        else
            X = X_clean;
        end
        % GWT parameters that need to be set separately
        GMRAopts.ManifoldDimension = dataset.d;
        GMRAopts.errorType = 'absolute';
        GMRAopts.precision  = 5e-3; % only for leaf nodes
        
    case 'MeyerStaircase-d-skew'
        dataset = struct();
        dataset.name = pExampleName;
        dataset.N = nSamples;
        if nargin>=3
            dataset.d = d;
        else
            dataset.d = 2;
        end
        dataset.D = 10000;
        dataset.MeyerStepWidth='variable';
        dataset.noiseLevel = 0/sqrt(dataset.D);
        
        % Generate data
        [X_clean,Labels,OutParams.GenerateDataSets] = GenerateDataSets( 'MeyerStaircase-d', ...
            struct('NumberOfPoints',dataset.N,'Dim',dataset.d,'width',dataset.MeyerStepWidth,...
            'EmbedDim',dataset.D,'NoiseType','Gaussian','NoiseParam',0) );
        
        % Add noise
        if dataset.noiseLevel>0,
            X = X_clean + dataset.noiseLevel*random('norm', 0,1, size(X_clean));
            GMRAopts.X_clean = X_clean;
        else
            X = X_clean;
        end
        % GWT parameters that need to be set separately
        GMRAopts.ManifoldDimension = dataset.d;
        GMRAopts.errorType = 'absolute';
        GMRAopts.precision  = 5e-3; % only for leaf nodes
        
        
    case 'Meyerstair'
        %% Meyerstair case
        w1      = 100;     % initial window width
        w2      = 5;       % final window width
        T1      = round(w1/5);  % translation step for the first point
        T2      = round(w2/5);  % translation step for the second to the last point
        X       = Generate_Meyer(nsample,w1,w2,T1,T2);
        GMRAopts.ManifoldDimension = 1;
        
    case 'D-Gaussian' % D-Gaussian
        %% data parameters
        dataset = struct();
        dataset.name = pExampleName;
        dataset.N = 100000;
        dataset.k = 0;
        if nargin>=3,
            dataset.D = d;
        end
        dataset.noiseLevel = 0.01/sqrt(dataset.D);
        
        % Generate data
        lFactor = 1/2;
        [X_clean,Labels,OutParams.GenerateDataSets] = GenerateDataSets( dataset.name, ...
            struct('NumberOfPoints',dataset.N,'Dim',dataset.D,'EmbedDim',dataset.D,'NoiseType','Gaussian','NoiseParam',0, ...
            'GaussianMean',[ones(5,1);lFactor^2*ones(10,1);lFactor^3*ones(20,1);lFactor^4*ones(40,1);zeros(dataset.D-75,1)]', ...
            'GaussianStdDev',0.2*[ones(5,1);lFactor^2*ones(10,1);lFactor^3*ones(20,1);lFactor^4*ones(40,1);zeros(dataset.D-75,1)]) ); % figure;plot(idct(X_clean(randi(size(X_clean,1),1),:)))
        
        % Add noise
        if dataset.noiseLevel>0,
            X = X_clean + dataset.noiseLevel*random('norm', 0,1, size(X_clean));
            GMRAopts.X_clean = X_clean;
        else
            X = X_clean;
        end
        
        %% GWT parameters that need to be set separately
        GMRAopts.ManifoldDimension = 0;
        GMRAopts.errorType = 'relative';
        GMRAopts.threshold0 = lFactor/2; % threshold for choosing pca dimension at each nonleaf node
        GMRAopts.precision  = lFactor/10; % only for leaf nodes
        
    case 'STL10-traing'
        
        % generate the dataset
        dataset = struct();
        dataset.digits = 0:9;
        dataset.projectionDimension = Inf;
        
        [X0,Labels,OutParams.GenerateDataSets] = GenerateDataSets( 'STL10traing', struct('AutotuneScales',false) );
        
        % image parameters
        imgOpts.imageData = true;
        imgOpts.imR = 96;
        imgOpts.imC = 96;
        imgOpts.Labels = Labels;
        
        if dataset.projectionDimension>0 && dataset.projectionDimension<imgOpts.imR*imgOpts.imC,                                % MM:TBD: to be changed, return projection, etc..as for MNIST
            imgOpts.X0 = X0;
            imgOpts.cm = mean(X0,2);
            X = X0 - repmat(imgOpts.cm,1,size(X0,2));
            [U,S,V] = randPCA(X, dataset.projectionDimension);
            X = S*V';
            imgOpts.U = U;
            imgOpts.isCompressed = true;
        else
            X = X0; clear X0;
            imgOpts.isCompressed = false;
        end
        
        
        % GWT parameters that need to be set separately
        %GWTopts.ManifoldDimension = 4; % if 0, then determine locally adaptive dimensions using the following fields:
        GMRAopts.threshold0 = 0.5; % threshold for choosing pca dimension at each nonleaf node
        GMRAopts.errorType = 'relative';
        GMRAopts.precision  = .1; % only for leaf nodes
    case 'Multiscale Patches'
        [imgOpts.FileName,imgOpts.PathName] = uigetfile({'*.jpg';'*.bmp';'*.png';'*.*'},'Pick an image to load');
        if imgOpts.FileName ~= 0
            GMRAopts.GenerateDataSets.FileName      = [imgOpts.PathName imgOpts.FileName];
            [X,Labels,OutParams.GenerateDataSets]    = GenerateDataSets('patchesfromimage', ...
                struct('FileName',GMRAopts.GenerateDataSets.FileName,'Labels',[],'AutotuneScales',false) );

            imgOpts.imageData   = true;
            imgOpts.imR         = OutParams.GenerateDataSets.GraphDiffOpts.imgSize(1);
            imgOpts.imC         = OutParams.GenerateDataSets.GraphDiffOpts.imgSize(2);
            imgOpts.Labels      = [];
            imgOpts.Coords      = OutParams.GenerateDataSets.Coords;
            
            GMRAopts.errorType  = 'relative';
            GMRAopts.precision  = .0001; % only for leaf nodes
            GMRAopts.ManifoldDimension = 0;
        end
    case 'Caltech101'
        dataset.projectionDimension = 100000;               
        %[X0,Labels,OutParams.GenerateDataSets] = GenerateDataSets( 'Caltech101', struct('Labels',{{'accordion','airplanes','brain','butterfly','camera','emu','hedgehog','lobster','pyramid','scissors'}},'SamplesPerClass',40,'NoiseParam',0,'AutotuneScales',false) );
        [X0,Labels,OutParams.GenerateDataSets] = GenerateDataSets( 'Caltech101', struct('Labels',{{'accordion','airplanes','hedgehog','scissors'}},'SamplesPerClass',100,'NoiseParam',0,'AutotuneScales',false) );
        ZeroIdx                 = find(sum(X0.^2,1)==0);
        X0(:,ZeroIdx)           = [];
        Labels(ZeroIdx)         = [];
        imgOpts.imageData       = true;
        imgOpts.imR             = OutParams.GenerateDataSets.GraphDiffOpts.imgSize(1);
        imgOpts.imC             = OutParams.GenerateDataSets.GraphDiffOpts.imgSize(2);
        imgOpts.Labels          = Labels;
        imgOpts.Coords          = OutParams.GenerateDataSets.Coords;
        
        if dataset.projectionDimension>0 && dataset.projectionDimension<imgOpts.imR*imgOpts.imC,
            if nargin<4 || ~isfield(GMRAopts,'Mean') || ~isfield(GMRAopts,'Proj') || isempty(GMRAopts.Mean) || isempty(GMRAopts.Proj)
                GMRAopts.Mean   = mean(X0,2);
                X               = bsxfun(@minus,X0,GMRAopts.Mean);
                [U,S,V]         = randPCA(X, dataset.projectionDimension);
                GMRAopts.Proj   = U';
            end
            X                   = GMRAopts.Proj*bsxfun(@minus,X0,GMRAopts.Mean);
            imgOpts.U           = GMRAopts.Proj;
            imgOpts.cm          = GMRAopts.Mean;
            imgOpts.isCompressed= true;
        else
            X = X0;
            GMRAopts.Proj = [];
            GMRAopts.Mean = [];
            imgOpts.isCompressed = false;
        end
        
    case 'Music'
        [X,Labels,OutParams.GenerateDataSets] = GenerateDataSets( 'Music', struct('AutotuneScales',false) );
    case 'Hyperspectral imaging'
        %load casi_2013.mat     xhsi = x2013_IEEE_GRSS_DF_Contest_CASI;
        load Indian_pines.mat;  xhsi = indian_pines;
        n1   = size(xhsi,1);
        n2   = size(xhsi,2);
        X    = zeros(size(xhsi,3),n1*n2);
        iter = 1;
        for k = 1:n1
            for j = 1:n2
                X(:,iter) = xhsi(k,j,:);
                iter      = iter+1;
            end
        end
    case '3D-surface'
        pointtype = 'dragon';
        switch pointtype
            case 'bunny'
               ptCloud = pcread('bun_zipper.ply');   %
               X       = (ptCloud.Location)';    %   
            case 'dragon'
                ptCloud = pcread('dragon.ply');  %pcread('xyzrgb_dragon.ply');   %
                X       = (ptCloud.Location)';    %   
            case 'teapot'
                ptCloud = pcread('teapot.ply');   %
                X       = (ptCloud.Location)';    % 
                X(3,:)      = X(3,:)-1;
            case 'armadillo'
                load armadillo01.mat
                X = (TR.Points)';
            case 'table'
                [vertex, ~] = read_off('b447.off');
                X           = vertex;
        end
        GMRAopts.ManifoldDimension = d;
    case 'whale'
        load mobysound_bow_hum_mfcc.mat
        X = H3';
        GMRAopts.ManifoldDimension = 0;
    case 'CIFAR10'
        LoadCifar10
        X = ConvertImagesToGrayScale(X,[32,32]); 
        X = abs(fft(X));
        Labels              = Y;
        X0                  = X;
        imgOpts.imageData   = true;
        imgOpts.imR         = 32;
        imgOpts.imC         = 32;
        imgOpts.Labels      = Y;
        imgOpts.isCompressed= false;
    case 'CIFAR100'
        %load('/Volumes/RAID/Data Sets/CIFAR/cifar-100-matlab/train.mat');
        load('../DataSets/CIFAR/cifar-100-matlab/train.mat');
        for k = size(data,1):-1:1
            tmp=reshape(data(k,:),[32,32,3]);
            tmp=sum(tmp,3);
            X0(:,k)=single(tmp(:))';
        end
        X = X0;
        imgOpts.imageData   = true;
        imgOpts.imR         = 32;
        imgOpts.imC         = 32;
        imgOpts.Labels      = [coarse_labels,fine_labels];
        GMRAopts.Proj = [];
        GMRAopts.Mean = [];
        imgOpts.isCompressed = false;
        
    case 'Faces BU3DFE'
        [X,ImInfo] = LoadFaces_BU3DFE;
        imgOpts.imageData   = true;
        imgOpts.imR         = ImInfo.imSz(1);
        imgOpts.imC         = ImInfo.imSz(1);
        imgOpts.Labels      = [];
        GMRAopts.Proj = [];
        GMRAopts.Mean = [];
        imgOpts.isCompressed = false;
        
    case 'PumaDyn'
       [X_train,Y_train,X_test,Y_test,Names] = pumadyn_CreateDataSets;
       
        lString = '\nWhich PumaDyn data set? \n';
        for k = 1:length(Names)
            lString = [lString '[' num2str(k) '] ' Names{k} '\n'];
        end
        while true
            k = input(lString);
            try
                if k>=1 && k<=length(X_train)
                    break;
                end
            end
        end
        
        X0     = [X_train{k},X_test{k}];
        Labels = [Y_train{k};Y_test{k}];        
        X      = X0;
        
    case 'HIGGS'
        load HIGGS                                                                                                              % This should be in the DataSets/HIGGS directory
        
        X0      = X;
        Labels  = [Y,Xderived'];
        
    otherwise % artificial data
        %% data parameters
        dataset = struct();
        dataset.name = pExampleName;
        if nargin == 1
            dataset.N = 100000;
        else
            dataset.N = nSamples;
        end
        if nargin>=3,
            dataset.k = d;
        end
        dataset.D = D;
        dataset.noiseLevel = 0/sqrt(dataset.D);
        
        % Generate data
        [X_clean,Labels,OutParams.GenerateDataSets] = GenerateDataSets( dataset.name, struct('NumberOfPoints',dataset.N,'Dim',dataset.k,'EmbedDim',dataset.D,'NoiseType','Uniform','NoiseParam',0,'AutotuneScales',false) );
        %[X_clean , Qrnd] = GenerateDataSets( dataset.name, struct('NumberOfPoints',dataset.N,'Dim',dataset.k,'EmbedDim',dataset.D,'NoiseType','Uniform','NoiseParam',0,'RandomProjectionDim',dataset.D) );
        %X_clean = X_clean';         % This is for 'Cosine' only.
        
        % Add noise
        if dataset.noiseLevel>0,
            X = X_clean + dataset.noiseLevel*random('norm', 0,1, size(X_clean)); % N by D data matrix
            GMRAopts.X_clean = X_clean;
        else
            X = X_clean;
        end
        
        %% GWT parameters that need to be set separately
        GMRAopts.ManifoldDimension = dataset.k;
        %GWTopts.threshold0=0.5;
        GMRAopts.errorType = 'absolute';
        GMRAopts.precision  = 1e-2; % only for leaf nodes
        
end

if nSamples<size(X,2)
    OutParams.AllX = X;
    perm           = randperm(size(X,2),nSamples);
    X              = X(:,perm);
    try X0         = X0(:,perm);  catch end
    if size(Labels,1)==1;  Labels = Labels'; end
    try Labels     = Labels(perm,:);  catch; end
end

if GMRAopts.verbose>0 
    fprintf('done. (%.3f sec)',toc);
end

% threshold for wavelet coefficients
GMRAopts.coeffs_threshold = 0; %GWTopts.precision/10;

%figure; do_plot_data(X,[],struct('view', 'pca'));
