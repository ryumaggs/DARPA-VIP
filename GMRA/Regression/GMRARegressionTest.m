function Yhat = GMRARegressionTest( gMRAReg, X_test, opts )

%
% function Yhat = GMRARegressionTest( gMRAReg, X_test, opts )
%
% IN:
%   gMRAReg     : GMRA regression structure as constructed by GMRARegression
%   X_test      : DxN matrix of data points
%   [opts]      : structure of options with the following fields:
%                   [partitions]    :
%                   [allscales]     : consider uniform partitions at all scales. Default: true, if opts.kappa not provided
%                   [kappa]         : k-vector, consider adaptive partitions with the threshold set to each of the values of kappa
%                   [Scores]        : Scores for X_test for polynomial regression. Default: empty, it is computed and returned in Yhat
%                   [ScalCoeffsExt] : Scaling function coefficients at each node of the GMRA. Default: empty, it is computed and returned in Yhat
%                   [Ypreds]        : predictions at each node. Default: empty, it is computed and returned in Yhat
%                   [XGWT]          : FGWT of X_test. Default: empty, it is computed and returned in Yhat
%
% OUT:
%   Yhat        : structure with the following fields:
%                   partition   : cell array of length equal to the number of polynomial degrees to try.
%                                 In the case of uniform partitions, partition{i} is a cell array whose l-th entry is the uniform partition at scale l
%                                 In the case of adaptive partitions, partition{i} is a cell array whose {l,r} entry is the adaptive partition corresponding to the
%                                   l-th threshold in the kappa{i,r} vector of thresholds for the r-th function
%                   Yhat        : cell array of length equal to the number of polynomial degrees to try
%
%
% (c) Mauro Maggioni, 2015
%

DEBUG = false;

if nargin<3,    opts = [];  end
if ~isfield(opts,'allscales') || ~opts.allscales
    if isfield(opts,'kappa')
        opts.allscales = false;
        if isempty(opts.kappa)
            MIN_KAPPA                   = 50;
            n                           = size(X_test,2);                                                                       % Construct a vector of reasonable values for kappa
            NKappas                     = 50;
            for polyidx = size(gMRAReg.Deltasq,1):-1:1                                                                          % Loop through the degrees of the local polynomials to fit
                for fcnidx = size(gMRAReg.Deltasq,3):-1:1                                                                       % Loop through the functions to be learned
                    for s = 1:10
                        tmp                     = 1/10*sqrt(gMRAReg.Deltasq(polyidx,:,fcnidx))/(sqrt(log(n)/n));
                        tmp2                    = tmp(tmp>max(tmp)*eps);
                        if isempty(tmp2), tmp2=tmp; end
                        Kappa{polyidx,fcnidx}                = 10.^unique(quantile(log10(tmp2),linspace(eps,1-eps,NKappas).^(1/2)));
                        if length(Kappa{polyidx,fcnidx})<MIN_KAPPA
                            NKappas = NKappas*2;
                        else
                            break
                        end
                    end
                    Kappa{polyidx,fcnidx}                    = sort(Kappa{polyidx,fcnidx},'descend');
                end
            end
            opts.kappa = Kappa;
        end
    end
else
    opts.allscales = true;
end
if ~isfield(opts,'Scores'),         opts.Scores         = [];    end
if ~isfield(opts,'ScalCoeffsExt'),  opts.ScalCoeffsExt  = [];    end
if ~isfield(opts,'Ypreds'),         opts.Ypreds         = [];    end
if ~isfield(opts,'XGWT') || isempty( opts.XGWT )                                                                                % Compute XGWT if not provided
    Yhat.XGWT           = FGWT( gMRAReg.gMRA,X_test );
else
    Yhat.XGWT           = opts.XGWT;
end
if isempty(opts.ScalCoeffsExt)
    Yhat.ScalCoeffsExt  = RearrangeScalingCoeffs( gMRAReg.gMRA, Yhat.XGWT );
else
    Yhat.ScalCoeffsExt  = opts.ScalCoeffsExt;
end
if isempty(opts.Scores)                                                                                                         % Compute regression scores if not provided
    Yhat.Scores         = ComputeRegressionScores( gMRAReg, X_test, Yhat.XGWT, Yhat.ScalCoeffsExt );
else
    Yhat.Scores         = opts.Scores;
end
if isempty(opts.Ypreds)
    for polyidx = 1:size(gMRAReg.Deltasq,1)
        [~,Yhat.Ypreds{polyidx}]  = PredictOnPartition( gMRAReg,polyidx,X_test,Yhat.XGWT,1:length(Yhat.XGWT.PointsInNet),Yhat.Scores, Yhat.ScalCoeffsExt );
    end
else
    Yhat.Ypreds         = opts.Ypreds;
end

% Construct regression on parts of the GMRA, processing top to bottom
[~,Yhat.nodeList]   = sort( gMRAReg.gMRA.Scales );

for polyidx = 1:size(gMRAReg.Deltasq,1)                                                                                         % Loop through degrees of polynomials
    if opts.allscales
        for j = min(gMRAReg.gMRA.Scales):max(gMRAReg.gMRA.Scales)
            Yhat.partition{polyidx}{j}   = get_partition_at_scale( gMRAReg.gMRA,j );
        end
    else
        Yhat.kappa  = opts.kappa;                                                                                               % Do thresholding depending on kappa
        partition   = cell(length(Yhat.kappa{polyidx,fcnidx}),size(gMRAReg.Deltasq,3));
        for fcnidx = 1:size(gMRAReg.Deltasq,3)                                                                                  % Loop through the functions
            for kappaidx = 1:length(Yhat.kappa{polyidx,fcnidx})                                                                 % Loop through the different values of kappa
                curdeltasq                         = squeeze(gMRAReg.Deltasq(polyidx,:,fcnidx));
                idxskeep                           = zeros(length(curdeltasq),1,'uint8');
                idxskeep(gMRAReg.gMRA.cp==0)       = 1;                                                        % Always keep at least the root
                idxskeep(curdeltasq>Yhat.kappa{polyidx,fcnidx}(kappaidx)*log(size(X_test,2))/size(X_test,2)) = 1;                
                [~,partition{kappaidx,fcnidx}]      = GetSmallestSubTreeWithIdxs( gMRAReg.gMRA.cp, idxskeep, gMRAReg.gMRA.pc );
            end
        end
        Yhat.partition{polyidx} = partition;
    end
end

Yhat.Yhat   = cell(size(gMRAReg.Deltasq,1),1);
for polyidx = 1:size(gMRAReg.Deltasq,1)                                                                                         % Loop through degrees of polynomials
    if opts.allscales
        for partitionidx = 1:size(Yhat.partition{polyidx},2)                                                                    % Loop through the partitions
            pred    = PredictOnPartition( gMRAReg, polyidx, X_test, Yhat.XGWT, Yhat.partition{polyidx}{partitionidx},Yhat.Scores, Yhat.ScalCoeffsExt, Yhat.Ypreds{polyidx} );
            if ~gMRAReg.RegressionOpts.isClassificationProblem
                Yhat.Yhat{polyidx}{partitionidx}   = pred;
            else
                [~,maxidx]                          = max(pred,[],2);
                Yhat.Yhat{polyidx}{partitionidx}   = gMRAReg.Labels(maxidx);
            end
        end
    else
        for fcnidx = 1:size(Yhat.partition{polyidx},2)                                                                          % Loop through the functions        
            for kappaidx = size(Yhat.partition{polyidx},1):-1:1                                                                 % Loop through partitions for each function
                pred = PredictOnPartition( gMRAReg, polyidx, X_test, Yhat.XGWT, Yhat.partition{polyidx}{kappaidx,fcnidx},Yhat.Scores, Yhat.ScalCoeffsExt, Yhat.Ypreds{polyidx}, fcnidx );
                if ~gMRAReg.RegressionOpts.isClassificationProblem
                    Yhat.Yhat{polyidx}{kappaidx}(:,fcnidx)   = pred;
                else
                    [~,maxidx]                                  = max(pred,[],2);               % MM: TBD!!!!!!!!!!!!!
                    Yhat.Yhat{polyidx}{kappaidx}(:,fcnidx)   = gMRAReg.Labels(maxidx);
                end
            end
        end
    end
end


return


function [Yhat,Yhat_local,idxs] = PredictOnPartition( gMRAReg, RegIdx, X, XGWT, Partition, Scores, ScalCoeffsExt, Ypreds, colidx )

DEBUG = false;

if nargin<9
    Yhat  = zeros(size(X,2),size(gMRAReg.Deltasq,3));
else
    Yhat  = zeros(size(X,2),length(colidx));
end

Yhat_local = cell(length(Partition),1);

idxs = cell(length(Partition),1);

if gMRAReg.RegressionOpts.noProj || isempty(ScalCoeffsExt)
    for k = 1:length(Partition)
        idxs{k}         = XGWT.PointsInNet{Partition(k)}';
        X_local         = X(:,idxs{k});
        if ~gMRAReg.RegressionOpts.noProj
            if ~isempty(gMRAReg.gMRA.ScalBasis{Partition(k)})
                X_local     = gMRAReg.gMRA.ScalBasis{Partition(k)}*bsxfun(@minus,X_local,gMRAReg.gMRA.Centers{Partition(k)});
                %ScalCoeffs = GWT_getScalCoeffs_atnode( gMRAReg.gMRA, XGWT, Partition(k) )';                                    % MM: this is actually slower
            end
            %figure;scatter3(X(1,idxs),X(2,idxs),X(3,idxs),20,log10(abs(Yhat(idxs)-fcn_norm(X_local))),'filled');colorbar;
        end
        Yhat_local{k}   = MultiPolyRegressEval( gMRAReg.Reg{RegIdx}{Partition(k)},X_local,Scores{RegIdx,Partition(k)} );
    end
else
    if nargin == 8
        Yhat_local = Ypreds(Partition);
    elseif nargin==9
        for k = 1:length(Partition)
            if ~isempty(Ypreds{Partition(k)})
                Yhat_local{k}       = Ypreds{Partition(k)}(:,colidx);
            end
        end
    else
        pointsInNet = XGWT.PointsInNet;
        for k = 1:length(Partition)
            idxs{k}             = pointsInNet{Partition(k)}';
            if ~isempty(ScalCoeffsExt{Partition(k)})
                Yhat_local{k}   = MultiPolyRegressEval( gMRAReg.Reg{RegIdx}{Partition(k)},ScalCoeffsExt{Partition(k)},Scores{RegIdx,Partition(k)} );
            else
                Yhat_local{k}   = MultiPolyRegressEval( gMRAReg.Reg{RegIdx}{Partition(k)},zeros(1,length(idxs{k})),Scores{RegIdx,Partition(k)} );
            end
        end
    end
end

for k = length(Partition):-1:1
    if ~isempty(Yhat_local{k,1})
        idxs{k}         = XGWT.PointsInNet{Partition(k)}';
        Yhat(idxs{k},:) = Yhat_local{k};
    end
end

return