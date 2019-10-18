function DataError = AGWT_ComputeApproxErr(gMRA,XGWT,X_test,XGWT_test,opts)
%
% function to compute approximation errors of uniform and adaptive GMRA
%
% DataError = AGWT_ComputeApproxErr(gMRA,XGWT,X_test,XGWT_test,opts)
%
% IN:
%   gMRA        : GMRA data structure
%   XGWT        : FGWT of training data
%   X_test      : test data
%   XGWT_test   : FGWT of test data
%   [opts]      : options:
%       [refinementnorm]: norm of refinement criterion: 
%                         '2': L^2 norm, 'inf': L^\infty norm
%       [ifscaled]: 1: scale-dependent threshold 0: scale-independent
%
% OUT:
%   DataError   : structure with error
%                 Let J be the maximal level of the tree
%       UniformRadii:              vector of length J whose jth entry is the average radius at level j 
%       UniformScales:             vector of length J whose jth entry is numeric scale  at level j 
%       UniformPartition:          uniform partition at level j = 1,..,J
%       UniformPartitionEntropy:   vector of length J whose jth entry is the entropy of the uniform partition at level j
%       UniformNcell:              vector of length J whose jth entry is the cardinality of the uniform partition at level j
%       UniformAbsolute:           vector of length J whose jth entry is the absolute error of GMRA at level j 
%       UniformRelative:           vector of length J whose jth entry is the relative error of GMRA at level j 
%       CenterAbsolute:            vector of length J whose jth entry is the absolute error of center approximation at level j 
%       CenterRelative:            vector of length J whose jth entry is the relative error of center approximation at level j 
%       AdaptiveKappa:             $\kappa$ in adaptive GMRA
%       AdaptiveAbsolute:          vector of absolute error of adaptive GMRA for the chosen $\kappa$ 
%       AdaptiveRelative:          vector of relative error of adaptive GMRA for the chosen $\kappa$      
%       AdaptivePartition:         adaptive partition for the chosen $\kappa$
%       AdaptivePartitionEntropy:  vector of entropy of the adaptive partition for the chosen $\kappa$ 
%   
%
%
% (c) Wenjing Liao, Mauro Maggioni


DEBUG = false;
if ~isfield(opts,'refinementnorm'), opts.refinementnorm = '2'; end         % default refinement norm: L2
if ~isfield(opts,'ifscaled'), opts.ifscaled = 1; end                       % default refinement: scale-dependent refinement

% Variables
J              = max(gMRA.Scales);  %max(find(min(XGWT_test.Cel_cpidx,[],1)));
MIN_KAPPA      = J;

GMRAErrOpts   = struct('norm',[2,inf],'relative',true,'quantiles',[],'type','GMRA'); %[0.25,0.5,0.75]);                                                   % Default option for approximation error calculation
CenterErrOpts = struct('norm',[2,inf],'relative',true,'quantiles',[],'type','Center'); %[0.25,0.5,0.75]);                                                   % Default option for approximation error calculation

% Uniform approximation errors
DataError.UniformRadii             = zeros(J,1);
DataError.UniformScales            = zeros(J,1);
DataError.CenterRadii              = zeros(J,1);
DataError.CenterScales             = zeros(J,1);
DataError.UniformScales            = zeros(J,1);
DataError.UniformPartitionEntropy  = zeros(J,1);
DataError.UniformTreeEntropy       = zeros(J,1);
DataError.UniformPartition         = cell(J,1);
DataError.CenterPartition          = cell(J,1);
%DataError.PercentBd                = zeros(J,1);

for j = 1:J
    DataError.UniformPartition{j}               = GetPartition(gMRA,XGWT,struct('type','Uniform','scale',j));                   % Get adaptive partition
    % compute the percentage of cells with size >= d at scale j
%     if gMRA.opts.ManifoldDimension > 0
%         CellSizes                               = gMRA.Sizes(DataError.UniformPartition{j}.nodes);
%         DataError.PercentBd(j)                  = length(find(CellSizes>=gMRA.opts.ManifoldDimension))/length(CellSizes); 
%     end
    % GMRA approximation
    [ApproxError_tmp, Quantiles_tmp]            = GetApproximationError(gMRA, X_test, XGWT_test, DataError.UniformPartition{j}.leafParentInPartition, GMRAErrOpts );
    DataError.UniformAbsolute(j,:)              = ApproxError_tmp.Absolute;
    DataError.UniformRelative(j,:)              = ApproxError_tmp.Relative;
    if ~isempty(Quantiles_tmp)
        DataError.UniformAbsoluteQuantiles(j,:,:)   = Quantiles_tmp.Absolute;
        DataError.UniformRelativeQuantiles(j,:,:)   = Quantiles_tmp.Relative;
    end
    DataError.UniformNCell(j)                   = length(DataError.UniformPartition{j}.nodes);
    DataError.UniformRadii(j)                   = mean(gMRA.Radii(DataError.UniformPartition{j}.nodes));
    [DataError.UniformPartitionEntropy(j),DataError.UniformTreeEntropy(j)]  = GetPartitionEntropy( gMRA, XGWT_test, DataError.UniformPartition{j} );                % Get entropy of partition
    % Center approximation
    DataError.CenterPartition{j}                = GetPartition(gMRA,XGWT,struct('type','Uniform','scale',j));                   % Get adaptive partition
    [ApproxError_tmp, Quantiles_tmp]            = GetApproximationError(gMRA, X_test, XGWT_test, DataError.CenterPartition{j}.leafParentInPartition, CenterErrOpts );
    DataError.CenterAbsolute(j,:)               = ApproxError_tmp.Absolute;
    DataError.CenterRelative(j,:)               = ApproxError_tmp.Relative;
    DataError.CenterRadii(j)                    = mean(gMRA.Radii(DataError.CenterPartition{j}.nodes));
end
if ~isfield(gMRA,'CoverTree') || ~isfield(gMRA.CoverTree,'theta'), gMRA.CoverTree.theta = 2;    end
DataError.UniformScales = log10(DataError.UniformRadii)/log10(gMRA.CoverTree.theta);
DataError.UniformScales(isinf(DataError.UniformScales)) = max(DataError.UniformScales);
DataError.CenterScales  = log10(DataError.CenterRadii)/log10(gMRA.CoverTree.theta);
DataError.CenterScales(isinf(DataError.CenterScales)) = max(DataError.CenterScales);
%[DataError,DataErrorQuantiles] = GetApproximationErrors(X_test, gMRA, XGWT, ErrOpts );                                         % DEBUG only

% Adaptive approximation error
n                           = size(XGWT.DeltaPoint,1);
NKappas                     = 32;
while true
    switch opts.refinementnorm
        case '2'
            if opts.ifscaled
                tmp                     = XGWT.DeltaCell./max(XGWT.CellRadii,min(min(XGWT.CellRadii(XGWT.CellRadii>0))))/(sqrt(log(n)/n));
            else
                tmp                     = XGWT.DeltaCell/(sqrt(log(n)/n)); 
            end
        case 'inf'
            if opts.ifscaled
                tmp                     = XGWT.DeltaCellInf./max(XGWT.CellRadii,min(min(XGWT.CellRadii(XGWT.CellRadii>0))))/(sqrt(log(n)/n)); 
            else
                tmp                     = XGWT.DeltaCellInf/(sqrt(log(n)/n));
            end
    end
    tmp                     = tmp(:);
    tmp                     = tmp(tmp>max(tmp)*10^(-6));
    Kappa                   = unique(quantile(tmp,linspace(0,1,NKappas)));
    %Kappa                   = union(Kappa,linspace(min(Kappa),max(Kappa),NKappas));
    if length(Kappa)<MIN_KAPPA,
        NKappas = NKappas*2;
    else
        break;
    end
end
Kappa                                = sort(Kappa,'descend');
LKappa                               = length(Kappa);
DataError.AdaptiveKappa              = Kappa';
DataError.AdaptiveAbsolute           = zeros(LKappa,2);
DataError.AdaptiveRelative           = zeros(LKappa,2);
DataError.AdaptiveNCell              = zeros(LKappa,1);
DataError.AdaptivePartitionEntropy   = zeros(LKappa,1);
DataError.AdaptiveTreeEntropy        = zeros(LKappa,1);

%% Truncate the tree and form adaptive partition
for j = 1 : LKappa
    kappa    = Kappa(j);
    DataError.AdaptivePartition{j}                  = GetPartition(gMRA,XGWT,struct('type','Adaptive','kappa',kappa,'disregardleaves',j~=LKappa,'refinementnorm',opts.refinementnorm,'ifscaled',opts.ifscaled));          % Get adaptive partition
    [ApproxError_tmp, Quantiles_tmp]                = GetApproximationError(gMRA, X_test, XGWT_test, DataError.AdaptivePartition{j}.leafParentInPartition, GMRAErrOpts );
    DataError.AdaptiveAbsolute(j,:)                 = ApproxError_tmp.Absolute;
    DataError.AdaptiveRelative(j,:)                 = ApproxError_tmp.Relative;
    if ~isempty(Quantiles_tmp)
        DataError.AdaptiveAbsoluteQuantiles(j,:,:)  = Quantiles_tmp.Absolute;
        DataError.AdaptiveRelativeQuantiles(j,:,:)  = Quantiles_tmp.Relative;
    end
    
    if (DEBUG)                                                                                                                  % DEBUG only
         DataAGWT(j)                                = FAGWT_TruncateTree(XGWT,kappa);                                                 
         if max(abs(DataAGWT(j).PartitionIndex-DataError.AdaptivePartition{j}.leafParentInPartition))>0, 
             warning('partitions are different!'), keyboard,    
         end
         X_Approx                                   = FAGWT_Approximation(gMRA,DataAGWT(j).PartitionIndex,XGWT_test);
         tdiff                                      = sum((X_test-X_Approx).^2,1);
         DataError.OldAdaptiveAbsError(j)           = sqrt(mean(tdiff));
         DataError.OldAdaptiveReError(j)            = sqrt(mean(tdiff./tx));    
         if abs(DataError.OldAdaptiveAbsError(j)-DataError.AdaptiveAbsolute(j,:))>100*eps
             warning('absolute errors are different!'), keyboard,    
         end
         if abs(DataError.OldAdaptiveReError(j)-DataError.AdaptiveRelative(j,:))>100*eps
             warning('relative errors are different!'), keyboard,    
         end
    end
    
    DataError.AdaptiveNCell(j)                  = length(DataError.AdaptivePartition{j}.nodes);
    [DataError.AdaptivePartitionEntropy(j) , DataError.AdaptiveTreeEntropy(j)] = GetPartitionEntropy( gMRA, XGWT_test, DataError.AdaptivePartition{j} );                % Get entropy of partition
end

return
