%% Script helping with the selection of various options

% (c) Copyright, Mauro Maggioni, Wenjing Liao

if ~exist('Params','var'),  Params = []; end

%% Pick a data set
ExampleNames  = {'MNIST_Digits','YaleB_Faces','croppedYaleB_Faces','NaturalImagePatches','ScienceNews', 'IntegralOperator','MeyerStaircase-d','MeyerStaircase-d-skew', ...
    'SwissRoll','SwissRollSkew','S-Manifold','Z-Manifold','Curved-Z-manifold','Oscillating2DWave','D-Ball','D-Sphere', 'D-Cube','D-FlatTorus','Cosine','Signal1D_1','D-Gaussian','STL10-traing', ...
    'Multiscale Patches','Caltech101','Music','CIFAR100','Faces BU3DFE','PumaDyn','HIGGS','Hyperspectral imaging','3D-surface','CIFAR10','Puma'};

fprintf('\n Examples:\n');
for k = 1:length(ExampleNames)
    fprintf('\n [%d] %s',k,ExampleNames{k});
end;
fprintf('\n\n');

if ~isfield(Params,'ExampleName')
    ExampleIdx          = input('Pick an example to run:           ');
    Params.ExampleName  = ExampleNames{ExampleIdx};
end

%% Choose number of points
if ~isfield(Params,'N')
    Params.N = input('Number of points [default=20000]:           ');
end
if isempty(Params.N) || isnan(Params.N), Params.N = 20000; end

%% Choose intrinsic dimension, when applicable
IfRealData  = [1 1 1 1 1    0 0 0 0 0  0 0 0 0 0   0 0 0 0 0    0 1 1 1 1  1 1 1 1 1   0 1 0];
if ~isfield(Params,'d')
    if exist('ExampleIdx','var') && ~IfRealData(ExampleIdx)
        if strcmp(Params.ExampleName,'Puma')
            njoints = input('Choose number of free joints [from 1 to 6, default=3]:           ');
            if floor(njoints)~=njoints || njoints<1 || njoints>6
                error('The number of joints should be an integer between 1 and 6.');
                return;
            end
            Params.d = 3*njoints-1;
            clear njoints;
        elseif ~isfield(Params,'d')
            Params.d = input('Pick a dimension [default=2 if manifold, appropriate d if not]:           ');
        end
        if isempty(Params.d), Params.d = 2; end
    else
        Params.d = 0;
    end
end

%% Choose noise level
if ~isfield(Params,'Sigma')
    Params.Sigma = input('Noise level on X (it will be divided by 1/sqrt(D)) [default=0]:          ');
end
if isempty(Params.Sigma), Params.Sigma = 0; end

%% Choose a GWT version
if ~exist('DIFFUSION_GEOMETRY_ONLY','var') && ~isfield(Params,'GWTversion')
    GWTversions  = {'Vanilla GWT','Orthogonal GWT','Pruning GWT'};
    methodLabels = [0 1 2];
    fprintf('\n Geometric Wavelets version:\n');
    for k = 1:length(GWTversions)
        fprintf('\n [%d] %s',methodLabels(k),GWTversions{k});
    end;
    fprintf('\n\n  ');
    
    if ~isfield(Params,'GWTversion')
        Params.GWTversion = input('Pick a version of the GWT to run [default=0]:           ');
    end
    if isempty(Params.GWTversion), Params.GWTversion = 0; end;
    
    %% Choose a partitioning method
    PartitionMethods  = {'METIS','Cover Tree'};
    methodLabels = [0 1];
    fprintf('\n Multiscale partitioning method:\n');
    for k = 1:length(PartitionMethods)
        fprintf('\n [%d] %s',methodLabels(k),PartitionMethods{k});
    end;
    fprintf('\n\n  ');
    
    if ~isfield(Params,'PartitionMethods')
        Params.PartitionMethods = input('Pick a multiscale partitioning method to use [default=1]:           ');
    end
    if isempty(Params.PartitionMethods), Params.PartitionMethods = 1; end
end

%% Choose if train & test
if ~isfield(Params,'TestPercent')
    Params.TestPercent = input('Size of test set in percentage [default=50]:           ');
end
if isempty(Params.TestPercent) 
    if ~exist('DIFFUSION_GEOMETRY_ONLY','var')
        Params.TestPercent = 50; 
    else
        Params.TestPercent = 0;
    end
end
