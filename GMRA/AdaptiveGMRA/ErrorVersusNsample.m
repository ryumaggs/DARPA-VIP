%% MSE versus #sample --- Compare Uniform and adaptive GMRA
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters
Times    = 5;
Nsample  = round(10.^(3:0.2:4+log10(5)));  %round(10.^(4:0.2:6+log10(5)));
LenN     = length(Nsample);
Kappa    = [0.1 0.5 1 2 4]; %[0.05 0.1 0.25 0.5 1 2];
LKappa   = length(Kappa);

%% Pick a data set
pExampleNames  = {'MNIST_Digits','YaleB_Faces','croppedYaleB_Faces','NaturalImagePatches','ScienceNews', 'IntegralOperator','MeyerStaircase-d','MeyerStaircase-d-skew', ...
    'SwissRoll','SwissRollSkew','S-Manifold','Z-Manifold','Curved-Z-manifold','Oscillating2DWave','D-Ball','D-Sphere', 'D-Cube','D-FlatTorus','Cosine','Signal1D_1','D-Gaussian','STL10-traing', ...
    'Multiscale Patches','Caltech101','Music','CIFAR100','Faces BU3DFE','PumaDyn','HIGGS','Hyperspectral imaging','Armadillo'};

fprintf('\n Examples:\n');
for k = 1:length(pExampleNames),
    fprintf('\n [%d] %s',k,pExampleNames{k});
end;
fprintf('\n\n  ');

pExampleIdx = input('Pick an example to run:           ');


%% Choose if train & test
pTestPercent = input('Size of test set in percentage [default=50]:           ');
if isempty(pTestPercent), pTestPercent = 50; end

%% Choose a GWT version
GWTversions  = {'Vanilla GWT','Orthogonal GWT','Pruning GWT'};
methodLabels = [0 1 2];
fprintf('\n Geometric Wavelets version:\n');
for k = 1:length(GWTversions),
    fprintf('\n [%d] %s',methodLabels(k),GWTversions{k});
end;
fprintf('\n\n  ');

pGWTversion = input('Pick a version of the GWT to run [default=0]:           ');
if isempty(pGWTversion), pGWTversion = 0; end;

%% Choose a partitioning method
PartitionMethods  = {'METIS','Cover Tree'};
methodLabels = [0 1];
fprintf('\n Multiscale partitioning method:\n');
for k = 1:length(PartitionMethods),
    fprintf('\n [%d] %s',methodLabels(k),PartitionMethods{k});
end; 
fprintf('\n\n  ');

pPartitionMethods = input('Pick a multiscale partitioning method to use [default=1]:           ');
if isempty(pPartitionMethods), pPartitionMethods = 1; end

%% Choose noise level
Sigma = input('\n Noise level on X (it will be divided by 1/sqrt(D)) [default=0]:          ');
if isempty(Sigma), Sigma = 0; end

%% Choose intrinsic dimension, when applicable
IfRealData  = [1 1 1 1 1    0 0 0 0 0  0 0 0 0 0   0 0 0 0 0    0 1 1 1 1  1 1 1 1 1   0];
if ~IfRealData(pExampleIdx)
    d = input('Pick a dimension [default=2 if manifold, appropriate d if not]:           ');
    if isempty(d), d = 2; end
else
    d = 0;
end


%% As and Bs
if pExampleIdx == 11
    AS = 2;   BS =2;
elseif pExampleIdx == 12
    AS = 1.5; BS =(3/2)*(d-2)/(d-3); 
elseif pExampleIdx == 13
    AS = 1.5;   BS = 1.5; 
else
    AS = 1;  BS =1;
end


IfSaveError = input('Save error data? 0: No; 1: Yes. [default = 1] \n');
IfSaveGMRA  = input('Save GMRA structure? 0: No; 1: Yes. [default = 0]\n');
if isempty(IfSaveError), IfSaveError=1; end
if isempty(IfSaveGMRA), IfSaveGMRA=0; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ScaleNames                       = {'K means','jopt',...
                                    'MedianSize>=d^2',...
                                    'jstar'};
LScale                           = length(ScaleNames);

DataErrorUniform0                = struct();
%
DataErrorUniform0.Allj           = zeros(Times,LenN,LScale);
DataErrorUniform0.Radiij         = zeros(Times,LenN,LScale);
DataErrorUniform0.Scalesj        = zeros(Times,LenN,LScale);
DataErrorUniform0.AbsErrAllj     = zeros(Times,LenN,LScale);
DataErrorUniform0.ReErrAllj      = zeros(Times,LenN,LScale);
%
DataErrorUniform0.Radii          = cell(Times,LenN);
DataErrorUniform0.Scales         = cell(Times,LenN);
DataErrorUniform0.PercentBd      = cell(Times,LenN);
DataErrorUniform0.Partition      = cell(Times,LenN);
DataErrorUniform0.AbsError       = cell(Times,LenN);
DataErrorUniform0.ReError        = cell(Times,LenN);
% Center
CenterErrorUniform0              = struct();
CenterErrorUniform0.Radii        = cell(Times,LenN);
CenterErrorUniform0.Scales       = cell(Times,LenN);
CenterErrorUniform0.AbsError     = cell(Times,LenN);
CenterErrorUniform0.ReError      = cell(Times,LenN);
CenterErrorUniform0.Partition    = cell(Times,LenN);

% Adaptive
DataErrorAdaptive0                = struct();
DataErrorAdaptive0.Kappa          = Kappa;
DataErrorAdaptive0.AbsError       = zeros(Times,LenN,LKappa);
DataErrorAdaptive0.ReError        = zeros(Times,LenN,LKappa);
DataErrorAdaptive0.Partition      = cell(Times,LenN,LKappa);
DataErrorAdaptive0.Par_Radii      = cell(Times,LenN,LKappa);
DataErrorAdaptive0.Par_Scales     = cell(Times,LenN,LKappa);
DataErrorAdaptive0.Par_Radii_mean = zeros(Times,LenN,LKappa);
DataErrorAdaptive0.Par_Radii_std  = zeros(Times,LenN,LKappa);

if IfSaveGMRA
    DataGMRA                     = cell(Times,LenN);
end

dirs = {'/gtmp/wjliao/','/ztmp/Wenjing/','/vtmp/Wenjing/'};
for itime = 1:length(dirs)
    if exist(dirs{itime},'dir'),
        dir = dirs{itime};
        break;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for isample = 1 : LenN
    nsample = Nsample(isample);
    pN      = nsample;
    for itime = 1 : Times
        fprintf('\n')
        display('====================================================')
        fprintf('Example = %6.0f d = %6.0f sigma = %6.2f \n',pExampleIdx,d,Sigma)
        fprintf('isample = %6.0f itime = %6.0f n = %6.0f \n',isample, itime, nsample)
        %% Generate data
        GenerateTrainAndTestData
        tx0_test                  = sum(X0_test.^2,1);
        tx0_test                  = tx0_test.*(tx0_test>0)+10^(-8).*(tx0_test==0);
        tx_test                   = sum(X_test.^2,1);
        tx_test                   = tx_test.*(tx_test>0)+10^(-8).*(tx_test==0);
        %% GMRA opts
        SetGMRAParameters                                                                                                               %% Set GMRA parameters
        GMRAopts.CoverTreeBuildOpts.numlevels = int32(100); 
        GMRAopts.CoverTreeTrimOpts.TrimValue  = d;                          % Go deep: useful for seing bias/variance trade-off
        GMRAopts.verbose                      = 0;
        GMRAopts.ComputeWavelets              = 0;        
        %% Construct GMRA
        gMRA       = GMRA(X_train, GMRAopts);
        %% FGWT
        XGWT       = FGWT(gMRA, X_train);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% clean test data
        XGWT0_test = FGWT(gMRA, X0_test);
        %% GMRA -- uniform
        J                                           = max(gMRA.Scales);
        %
        DataErrorUniform0.theta(itime,isample)      = theta;
        DataErrorUniform0.PercentBd{itime,isample}  = zeros(J,1);       
        % K means error
        CenterErrorUniform0.Radii{itime,isample}    = zeros(1,J);
        CenterErrorUniform0.Scales{itime,isample}   = zeros(1,J);
        CenterErrorUniform0.AbsError{itime,isample} = zeros(J,1);
        CenterErrorUniform0.ReError{itime,isample}  = zeros(J,1); 
        CenterErrorUniform0.Partition{itime,isample}= cell(J,1);
        % Maximum partition size
        Max_part_size                               = zeros(1,J);
        %% Partition on the master tree whose leaf size = 1
        for j = 1:J
            % GMRA
            Partitionj                                      = GetPartition(gMRA,XGWT,struct('type','Uniform','scale',j,'minleafsize',1));  % partition at scale j
            CenterErrorUniform0.Partition{itime,isample}{j} = Partitionj;
            % K means            
            X0_test_Approx                                  = GetApproximationOnPartition(gMRA,Partitionj.leafParentInPartition,XGWT0_test,struct('type','Center')); % approximation at scale j
            % error
            tdiff                                           = sum((X0_test-X0_test_Approx).^2,1);
            CenterErrorUniform0.AbsError{itime,isample}(j)  = sqrt(mean(tdiff));
            CenterErrorUniform0.ReError{itime,isample}(j)   = sqrt(mean(tdiff./tx0_test));
            % Radii
            CenterErrorUniform0.Radii{itime,isample}(j)     = mean(gMRA.Radii(Partitionj.nodes)); 
            % Cell size at scale j
            CellSizes                                       = gMRA.Sizes(Partitionj.nodes);
            % percentage of points whose size >= d
            DataErrorUniform0.PercentBd{itime,isample}(j)   = length(find(CellSizes>=d))/length(CellSizes);          
            % minimal cell size at scale j 
            Max_part_size(j)                                = max(CellSizes); 
        end
        CenterErrorUniform0.Scales{itime,isample}           = log10(CenterErrorUniform0.Radii{itime,isample})/log10(theta);
        %% Partition on the master tree whose leaf size >= d
        Jd = find(Max_part_size>=d,1,'last');  % the maximum scale at which there exists a cell >= d
        % GMRA error on the master tree whose leaf size = 1
        DataErrorUniform0.Radii{itime,isample}      = zeros(1,Jd);
        DataErrorUniform0.Partition{itime,isample}  = cell(Jd,1); 
        DataErrorUniform0.AbsError{itime,isample}   = zeros(Jd,1);
        DataErrorUniform0.ReError{itime,isample}    = zeros(Jd,1);       
        Median_part_size                            = zeros(1,Jd); % Median partition size
        for j = 1 : Jd
            % GMRA
            Partitionj                                         = GetPartition(gMRA,XGWT,struct('type','Uniform','scale',j,'sigma',Sigma,'minleafsize',d));  % partition at scale j
            DataErrorUniform0.Partition{itime,isample}{j}      = Partitionj;
            X0_test_Approx                                     = GetApproximationOnPartition(gMRA,Partitionj.leafParentInPartition,XGWT0_test); % approximation at scale 
            % Error
            tdiff                                              = sum((X0_test-X0_test_Approx).^2,1);
            DataErrorUniform0.AbsError{itime,isample}(j)       = sqrt(mean(tdiff));
            DataErrorUniform0.ReError{itime,isample}(j)        = sqrt(mean(tdiff./tx0_test));
            % Radii
            DataErrorUniform0.Radii{itime,isample}(j)          = mean(gMRA.Radii(Partitionj.nodes)); 
            % Cell size at scale j
            CellSizes                                          = cellfun(@length,gMRA.PointsInNet(Partitionj.nodes));
            % Median cell size at scale j 
            Median_part_size(j)                                = median(CellSizes); 
        end
        DataErrorUniform0.Scales{itime,isample}                = log10(DataErrorUniform0.Radii{itime,isample})/log10(theta);
        % Record error at certain scales
                                    %{'K means','jopt',...
                                    %'jopt&MedianSize>=d^2',...
                                    %'jstar'};
        fprintf('\n clean test samples  Jd = %6.0f error at Jd =%6.6f\n',Jd,DataErrorUniform0.AbsError{itime,isample}(Jd));
        for iscale = 1:LScale
            scalename = ScaleNames{iscale};
            if strcmp(scalename,'K means')
                jchosen                                            = J;
                DataErrorUniform0.Allj(itime,isample,iscale)       = jchosen;
                DataErrorUniform0.Radiij(itime,isample,iscale)     = CenterErrorUniform0.Radii{itime,isample}(jchosen);
                DataErrorUniform0.Scalesj(itime,isample,iscale)    = CenterErrorUniform0.Scales{itime,isample}(jchosen);
                DataErrorUniform0.AbsErrAllj(itime,isample,iscale) = CenterErrorUniform0.AbsError{itime,isample}(jchosen);
                DataErrorUniform0.ReErrAllj(itime,isample,iscale)  = CenterErrorUniform0.ReError{itime,isample}(jchosen);
                fprintf([scalename '\n%6.0f  Abs Err %6.6f\n'],jchosen,CenterErrorUniform0.AbsError{itime,isample}(jchosen))
            elseif strcmp(scalename,'jopt')
                [~ ,jchosen]                                       = min(DataErrorUniform0.AbsError{itime,isample});
                DataErrorUniform0.Allj(itime,isample,iscale)       = jchosen;
                DataErrorUniform0.Radiij(itime,isample,iscale)     = DataErrorUniform0.Radii{itime,isample}(jchosen);
                DataErrorUniform0.Scalesj(itime,isample,iscale)    = DataErrorUniform0.Scales{itime,isample}(jchosen);
                DataErrorUniform0.AbsErrAllj(itime,isample,iscale) = DataErrorUniform0.AbsError{itime,isample}(jchosen);
                DataErrorUniform0.ReErrAllj(itime,isample,iscale)  = DataErrorUniform0.ReError{itime,isample}(jchosen);
                fprintf([scalename '\n%6.0f  Abs Err %6.6f\n'],jchosen,DataErrorUniform0.AbsError{itime,isample}(jchosen))
            elseif strcmp(scalename,'MedianSize>=d^2')
                jchosen                                            = find(Median_part_size>=d^2,1,'last');
                DataErrorUniform0.Allj(itime,isample,iscale)       = jchosen;
                DataErrorUniform0.Radiij(itime,isample,iscale)     = DataErrorUniform0.Radii{itime,isample}(jchosen);
                DataErrorUniform0.Scalesj(itime,isample,iscale)    = DataErrorUniform0.Scales{itime,isample}(jchosen);
                DataErrorUniform0.AbsErrAllj(itime,isample,iscale) = DataErrorUniform0.AbsError{itime,isample}(jchosen);
                DataErrorUniform0.ReErrAllj(itime,isample,iscale)  = DataErrorUniform0.ReError{itime,isample}(jchosen);
                fprintf([scalename '\n%6.0f  Abs Err %6.6f\n'],jchosen,DataErrorUniform0.AbsError{itime,isample}(jchosen))
            elseif strcmp(scalename,'jstar')
                [~ , jchosen]   = min(abs(DataErrorUniform0.Radii{itime,isample}- (log(nsample)/nsample)^(1/(2*AS+d-2))));   %max([2 , find(min_part_size>d^2,1,'last')]);
                jchosen         = min([jchosen Jd find((DataErrorUniform0.Radii{itime,isample}).^2>Sigma^2,1,'last')]);
                DataErrorUniform0.Allj(itime,isample,iscale)       = jchosen;
                DataErrorUniform0.Radiij(itime,isample,iscale)     = DataErrorUniform0.Radii{itime,isample}(jchosen);
                DataErrorUniform0.Scalesj(itime,isample,iscale)    = DataErrorUniform0.Scales{itime,isample}(jchosen);
                DataErrorUniform0.AbsErrAllj(itime,isample,iscale) = DataErrorUniform0.AbsError{itime,isample}(jchosen);
                DataErrorUniform0.ReErrAllj(itime,isample,iscale)  = DataErrorUniform0.ReError{itime,isample}(jchosen);
                fprintf([scalename '\n%6.0f  Abs Err %6.6f\n'],jchosen,DataErrorUniform0.AbsError{itime,isample}(jchosen))
            end
        end
        %% Adaptive GMRA
        for ikappa = 1:LKappa
            kappa = Kappa(ikappa);
            % Get adaptive partition
            AdaptivePartition                                      = GetPartition(gMRA,XGWT,struct('type','Adaptive','kappa',kappa,'sigma',Sigma,'minleafsize',d));          % Get adaptive partition
            % Approximation
            X0_test_Approx                                         = GetApproximationOnPartition(gMRA,AdaptivePartition.leafParentInPartition,XGWT0_test);
            % Error
            tdiff                                                  = sum((X0_test-X0_test_Approx).^2,1);                 
            DataErrorAdaptive0.AbsError(itime,isample,ikappa)      = sqrt(mean(tdiff));
            DataErrorAdaptive0.ReError(itime,isample,ikappa)       = sqrt(mean(tdiff./tx0_test));
            DataErrorAdaptive0.Partition{itime,isample,ikappa}     = AdaptivePartition;
            DataErrorAdaptive0.Par_Radii{itime,isample,ikappa}     = gMRA.Radii(AdaptivePartition.nodes);
            DataErrorAdaptive0.Par_Radii_mean(itime,isample,ikappa)= mean(DataErrorAdaptive0.Par_Radii{itime,isample,ikappa});
            DataErrorAdaptive0.Par_Radii_std(itime,isample,ikappa) = std(DataErrorAdaptive0.Par_Radii{itime,isample,ikappa});
            DataErrorAdaptive0.Par_Scales{itime,isample,ikappa}    = log10(DataErrorAdaptive0.Par_Radii{itime,isample,ikappa})/log10(theta);
        end
        %% save gMRA structure
        if IfSaveGMRA
            DataGMRA{itime,isample}                     = struct();
            gMRA                                        = rmfield(gMRA,{'Timing','sorted_idxs','cp_nIdxs','IniLabelsPtIdxs',...
                                                            'PointsInNet','WavSingVals','cp_orig','WavBases',...
                                                            'WavConsts','ScalCenter','ScalBasisChange','ScalWavBasisChange','ScalCenterConsts',...
                                                            'ScalBasisCoarserCenter','WavConstsCoeffs','DictCosts'});
            gMRA                                        = rmfield(gMRA,{'opts','CoverTree','isaleaf','Centers','ScalBasis',...
                                                             'Sigmas','TreeOnLeaves'});            
            DataGMRA{itime,isample}.gMRA                = gMRA;  
%             DataGMRA{itime,isample}.X_train             = X_train;
%             DataGMRA{itime,isample}.X0_test             = X0_test;
%             DataGMRA{itime,isample}.X_test              = X_test;
            XGWT                                        = rmfield(XGWT,{'Timing','dists2NearestLeafNodes','CelWavCoeffs',...
                                                               'CelScalCoeffs','PointsInNet','CelTangCoeffs','DeltaMatrix','CoeffsCosts',...
                                                               'CellMeasure','DeltaPoint'});
            XGWT0_test                                  = rmfield(XGWT0_test,{'Timing','dists2NearestLeafNodes','CelWavCoeffs',...
                                                               'CelScalCoeffs','PointsInNet','CelTangCoeffs','DeltaMatrix','CoeffsCosts',...
                                                               'CellMeasure','DeltaPoint'});
%           XGWT_test                                   = rmfield(XGWT_test,{'Timing','dists2NearestLeafNodes','CelWavCoeffs','CelTangCoeffs','DeltaMatrix','CoeffsCosts'});
            DataGMRA{itime,isample}.XGWT                = XGWT;
            DataGMRA{itime,isample}.XGWT0_test          = XGWT0_test;
%           DataGMRA{itime,isample}.XGWT_test           = XGWT_test;
        end
        %% save approximation error
        if IfSaveError
            save([dir 'Error8' pExampleNames{pExampleIdx} 'Dim' num2str(d) 'Sigma' num2str(Sigma) '.mat'],...
            'pExampleIdx','pExampleNames','d','Sigma','Nsample','DataErrorUniform0',...
            'CenterErrorUniform0','Kappa','DataErrorAdaptive0',...,
            'ScaleNames','IfRealData','-v7.3')
            fprintf('Error was saved in %s\n',dir)
        end
    end
end

DataErrorUniform0.Radiij_mean      =  mean(DataErrorUniform0.Radiij,1);
DataErrorAdaptive0.Par_Radii_mean2 =  mean(DataErrorAdaptive0.Par_Radii_mean,1);

DataErrorUniform0.AbsErrAllj_mean  =  mean(DataErrorUniform0.AbsErrAllj,1);
DataErrorAdaptive0.AbsError_mean   =  mean(DataErrorAdaptive0.AbsError,1);

%% save approximation error
if IfSaveError
    save([dir 'Error8' pExampleNames{pExampleIdx} 'Dim' num2str(d) 'Sigma' num2str(Sigma) '.mat'],...
            'pExampleIdx','pExampleNames','d','Sigma','Nsample','DataErrorUniform0',...
            'CenterErrorUniform0','Kappa','DataErrorAdaptive0','ScaleNames','IfRealData','-v7.3')
    fprintf('Error was saved in %s\n',dir)
end


 
%% save GMRA structure
if IfSaveGMRA
    save([dir 'GMRA8' pExampleNames{pExampleIdx} 'Dim' num2str(d) 'Sigma' num2str(Sigma) '.mat'],...
    'DataGMRA','GMRAopts','-v7.3')
end

