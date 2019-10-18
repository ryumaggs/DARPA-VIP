%%
% This file runs several examples of diffusion maps

clear Params

%% Selection of data set, and data generation/loading
DIFFUSION_GEOMETRY_ONLY = true;
SelectDataSetAndOptions                                                                                                         %% Select data set and other options
GenerateTrainAndTestData                                                                                                        %% Generate train and test samples

%% Diffusion Geometry Construction
fprintf('\n\nConstructing graph and diffusion map...\n');
GraphDiffOpts                   = OutParams.GenerateDataSets.GraphDiffOpts;
Data.G                          = GraphDiffusion(X, 0, GraphDiffOpts);                                                
fprintf('done.\n');

simpleVizForDiffusionGeometry

return

%
% This file runs several examples of diffusion maps
% If variables are not cleared before running this script and pExampleIdx exists, it will run the corresponding example without querying the user
%

clear pExampleIdx;close all;
clc

global X Data;

pExampleNames   = { '2-d sphere with little noise', ...
    '2-d sphere with large noise', ...
    'Sphere and segment', ...
    'Spiral and plane', ...
    'Meyerstaircase', ...
    'D-Cube with low sampling rate', ...
    'S-Manifold', ...
    '10-d cube with high noise', ...
    '9-d sphere with high noise', ...
    'Two lines and a plane', ...
    'Isomap faces', ...
    'CBCL faces, I', ...
    'CBCL faces, II', ...
    'Science News articles', ...
    'MNIST digits', ...
    'Patches from Lena Image'};

fprintf('\n\n Select example to run:\n');
for k = 1:length(pExampleNames),
    fprintf('\n [%d] %s',k,pExampleNames{k});
end;
fprintf('\n\n  ');

while true,
    if (~exist('pExampleIdx') || isempty(pExampleIdx) || pExampleIdx==0),
        try
            pExampleIdx = input('');
            pExampleIdx = str2num(pExampleIdx);
        catch
        end;
    end;
    if (pExampleIdx>=1) && (pExampleIdx<=length(pExampleNames)),
        break;
    else
        fprintf('\n %d is not a valid Example. Please select a valid Example above.',pExampleIdx);
        pExampleIdx=0;
    end;
end;

%% Set parameters for constructing the graph
EstDimOpts = struct('NumberOfTrials',15,'verbose',0,'MAXDIM',100,'MAXAMBDIM',100,'Ptwise',false,'NetsOpts',[],'UseSmoothedS',false, 'EnlargeScales',true );

Npts = 10000;

%% Set parameters for data generation and generates the data set
switch pExampleIdx
    case 1
        XName = 'D-Sphere'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',Npts,'Dim',2,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.1/sqrt(100));
        IsManifold = true;
    case 2
        XName = 'D-Sphere'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',Npts,'Dim',2,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.1);
        IsManifold = true;
    case 3
        XName = 'SphereAndLine'; XNickName = 'Sphere and Line';
        XOpts = struct('Dim','2.5','NoiseType','Gaussian','NoiseParam',0.00);
        EstDimOpts.Ptwise = true;
        IsManifold = false;
    case 4
        XName = 'SpiralAndPlane'; XNickName = 'Spiral and Plane';
        XOpts = struct('NumberOfPoints',1100,'Dim','1.5');
        EstDimOpts.Ptwise = true;
        Labels=[ones(300,1);2*ones(800,1)];
        IsManifold = false;
    case 5
        XName = 'MeyerStaircase'; XNickName = 'Z';
        XOpts = struct('NumberOfPoints',Npts,'Dim',1000,'MeyerStepWidth',20,'EmbedDim',1000,'NoiseType','Gaussian','NoiseParam',0.05/sqrt(1000));
        IsManifold = true;
    case 6
        XName = 'D-Cube'; XNickName = 'Q';
        XOpts = struct('NumberOfPoints',round(Npts/4),'Dim',6,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.01/sqrt(100));
        IsManifold = true;        
    case 7
        XName = 'S-Manifold'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',Npts,'Dim',2,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.01);
        IsManifold = true;
    case 8
        XName = 'D-Cube'; XNickName = 'Q';
        XOpts = struct('NumberOfPoints',Npts,'Dim',10,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.1);
        IsManifold = true;
    case 9
        XName = 'D-Sphere'; XNickName = 'S';
        XOpts = struct('NumberOfPoints',Npts,'Dim',9,'EmbedDim',100,'NoiseType','Gaussian','NoiseParam',0.1);
        EstDimOpts.EnlargeScales = true;
        IsManifold = true;
    case 10
        XName = 'TwoLinesAndAPlane'; XNickName = 'Two Lines and a Plane';
        XOpts = struct('NumberOfPoints',2*round(Npts/2),'Dim','1.5');
        EstDimOpts.Ptwise = true;
        Labels= [ones(round(Npts/2),1);2*ones(round(Npts/2),1)];
        IsManifold = false;
    case 11
        XName = 'Isomapfaces'; XNickName = 'Isomap Faces';
        XOpts = struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
        IsManifold = false;
    case 12
        XName = 'CBCLFaces1'; XNickName = 'CBCLFaces1 Faces';
        XOpts = struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
        IsManifold = false;
    case 13
        XName = 'CBCLFaces2'; XNickName = 'CBCLFaces2 Faces';
        XOpts = struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
        IsManifold = false;
    case 14
        XName = 'ScienceNews'; XNickName = 'ScienceNews Articles';
        XOpts = struct('BMarkUnitBall',true); %struct('SVDProject',200);
        EstDimOpts.Ptwise = true;
        IsManifold = false;
    case 15
        XName = 'BMark_MNIST'; XNickName = 'MNIST Digits';
        XOpts = struct('NumberOfPoints',200,'MnistOpts',struct('Sampling', 'RandN', 'QueryDigits',0:9, 'ReturnForm', 'vector'));
        IsManifold = false;
    case 16
        XName = 'ImagePatches'; XNickName = 'Lena Patches';
        XOpts = struct('ImageFileName','Lena.jpg','PatchSize',16,'DownsampleImage',8);
        IsManifold = false;
end

%% Generate the data set
fprintf('\nGenerating %s data...', XName);
XOpts.AutotuneScales    = false;
[X,Labels,OutParams]    = GenerateDataSets( XName, XOpts);
[Y]                     = GenerateDataSets( XName, XOpts);                                                                      % Create data where to test out-of-sample extension of diffusion map
X = 1*X;                                                                                                                        % Scale data if desired
fprintf('done.');

%% Compute Graph Diffusion and EigenFunctions
fprintf('\n\nConstructing graph and diffusion map...\n');
GraphDiffOpts                       = OutParams.GraphDiffOpts;
%GraphDiffOpts.DontReturnDistInfo    = 0;                                                                                        % These are for debugging purposes
if IsManifold, 
    GraphDiffOpts.Normalization     = 'beltrami'; 
    Data.G                          = GraphDiffusion(X, 0.5, GraphDiffOpts);                                                
else
    Data.G                          = GraphDiffusion(X, 0, GraphDiffOpts);                                                
end

fprintf('done.\n');

%% Simple visualization
% Original data and diffusion embeddings
scrsz = get(groot,'ScreenSize');
figure('Position',[scrsz(3)*1/8,scrsz(4)*1/8,scrsz(3)*3/4,scrsz(4)*3/4]);
subplot(1,3,1);plot3(X(1,:),X(2,:),X(3,:),'.');
xlabel('x_1','Interpreter','tex');ylabel('x_2','Interpreter','tex');zlabel('x_3','Interpreter','tex');
title('Original data, first 3 coords');
subplot(1,3,2);plot3(Data.G.EigenVecs(:,2),Data.G.EigenVecs(:,3),Data.G.EigenVecs(:,4),'.');
xlabel('\phi_2','Interpreter','tex');ylabel('\phi_3','Interpreter','tex');zlabel('\phi_4','Interpreter','tex');
title('Diffusion embedding (2,3,4)');
subplot(1,3,3);plot3(Data.G.EigenVecs(:,5),Data.G.EigenVecs(:,6),Data.G.EigenVecs(:,7),'.');
xlabel('\phi_5','Interpreter','tex');ylabel('\phi_6','Interpreter','tex');zlabel('\phi_7','Interpreter','tex');
title('Diffusion embedding (5,6,7)');

% Original data colored by values of eigenfunctions
figure('Position',[scrsz(3)*1/8,scrsz(4)*1/8,scrsz(3)*3/4,scrsz(4)*3/4]);
for k = 1:15
    subplot(3,5,k)
    scatter3(X(1,:),X(2,:),X(3,:),20,Data.G.EigenVecs(:,k),'filled');
    title(sprintf('Eigenfunction %d on the data',k));
end


fprintf('\n');
