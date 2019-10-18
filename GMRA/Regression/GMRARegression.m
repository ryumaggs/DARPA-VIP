function gMRAReg = GMRARegression( X, Y, opts )

%
% function gMRAReg = GMRARegression( X, Y, opts )
%
% IN:
%   X       : D by N matrix of N points in R^D
%   Y       : N by nY matrix of function values for nY different functions on X
%   [opts]  : structure of options with the following fields:
%       [degree] : vector of degrees of polynomials for local fitting. It will do a fit for every degree specified. Default: [1].
%       [noProj] : do not project on GMRA scaling subspaces. Default: false.
%       [gMRAReg]: data structure as returned by this function. If provided it will assume it was constructed on the same X
%                   and it will reuse the gMRA structure in it, as well as the local polynomial scores in order to speed up
%                   the calculation of the regression function. Default: [].
%       [isClassificationProblem] : construct classifier instead of regression. Default: false
% OUT:
%   gMRAReg : structure with teh following fields:
%           gMRA    : GMRA data structure for X. If opts.gGMRAReg is provided, it uses opts.gGMRAReg.gMRA
%           XGWT    : FGWT of the data X onto gMRA. If opts.gGMRAReg is provided, it uses opts.gGMRAReg.XGWT
%           Reg     : a cell array of length(opts.degree), with the i-th entry, corresponding to local polyonmials of degree opts.degree(i)x
%                       being a cell array of length equal to (#nodes in gMRA), each cell of which contains the local regression information
%                       (regression coefficients, scores, etc..., as returned by MultiPolyRegress).
%           Deltasq : a tensor of size length(opts.degree) x (#nodes in gMRA) x (size(Y,2)) with the squared error differences in approximation errors
%                       between each node and its parent
%           RegressionOpts : opts used
%           Labels  : class labels if opts.isClassificationProblem==true.
%
%
% (c) Mauro Maggioni, 2015
%

% Set options
if nargin<3,    opts = [];                                                          end
if ~isfield(opts,'degree'), opts.degree = 1;                                        end                                         % Can be a vector
if ~isfield(opts,'noProj'), opts.noProj = false;                                    end
if ~isfield(opts,'isClassificationProblem'), opts.isClassificationProblem = zeros(1,size(Y,2));  end
if isfield(opts,'gMRAReg'), gMRAReg = opts.gMRAReg;                                 end

warning off;

gMRAReg.Timings.overall = cputime;
gMRAReg.Timings.overall_tic = tic;

if opts.isClassificationProblem
    if size(Y,2)>1, warning('Run only one classification problem at a time'); Y=Y(:,1); end
    gMRAReg.Labels  = unique(Y);
    newY            = zeros(size(Y,1),length(gMRAReg.Labels));
    for k = 1:size(newY,2)
        newY(Y==gMRAReg.Labels(k),k) = 1;
    end
    Y = newY; clear newY;
end

if ~isfield(gMRAReg,'gMRA')
    % Construct GMRA on X
    if ~isfield(opts,'gMRA')
        gMRAReg.gMRA = GMRA( X,opts.GMRAopts );
    else
        gMRAReg.gMRA = opts.gMRA;
    end
    
    % Compute FGWT     
    gMRAReg.XGWT = FGWT( gMRAReg.gMRA,X );
    
    % Estimate function by local polynomials on the C_{j,k}'s of the GMRA
    for i = 1:length(opts.degree)
        gMRAReg.Reg{i}           = cell(length(gMRAReg.XGWT.PointsInNet),1);
    end
end

for i = 1:length(opts.degree)
    for k = length(gMRAReg.XGWT.PointsInNet):-1:1
        if isempty(gMRAReg.XGWT.PointsInNet{k}), continue; end;
        idxs                = gMRAReg.XGWT.PointsInNet{k};
        X_local             = X(:,idxs);
        if ~opts.noProj
            if ~isempty(gMRAReg.gMRA.ScalBasis{k})
                X_local         = gMRAReg.gMRA.ScalBasis{k}*bsxfun(@minus,X_local,gMRAReg.gMRA.Centers{k});
            else
                X_local         = zeros(1,length(idxs));
            end
        end
        Y_local             = Y(idxs,:);
        gMRAReg.Reg{i}{k}   = MultiPolyRegress( X_local,Y_local,opts.degree(i), 'reg',gMRAReg.Reg{i}{k} );                      % Fit local polynomial, possibly multiple degrees
    end
end

maxScale = max(gMRAReg.gMRA.Scales);
for i = 1:length(opts.degree)
    for j = maxScale:-1:1
        partitionj = get_partition_at_scale( gMRAReg.gMRA,j );
        for k = 1:length(partitionj)
            if ~isempty(gMRAReg.Reg{i}{partitionj(k)}.Coefficients)
                Yj{i}{j}(gMRAReg.XGWT.PointsInNet{partitionj(k)},:) = gMRAReg.Reg{i}{partitionj(k)}.yhat;
            end
        end
    end
end
gMRAReg.Yj = Yj;

% Compute the \Delta^2_{j,k}'s
[~,nodeList]            = sort(gMRAReg.gMRA.Scales);                                                                            % List nodes and their parents in top-bottom (coarsest to finest) fashion
parentList              = gMRAReg.gMRA.cp(nodeList);
parentList              = parentList(parentList>0);
gMRAReg.Deltasq         = zeros(length(opts.degree),length(nodeList),size(Y,2));
for k = 1:length(parentList)
    if parentList(k)>0
        for i = 1:length(opts.degree)
            gMRAReg.Deltasq(i,parentList(k),:) = gMRAReg.gMRA.Radii(parentList(k))*sum((gMRAReg.Yj{i}{gMRAReg.gMRA.Scales(parentList(k))}(gMRAReg.XGWT.PointsInNet{parentList(k)},:) - ...  % Multiply by gMRAReg.gMRA.Radii(parentList(k)) for Lip functions?
                gMRAReg.Yj{i}{gMRAReg.gMRA.Scales(parentList(k))+1}(gMRAReg.XGWT.PointsInNet{parentList(k)},:)).^2);
            %gMRAReg.Reg{i}{parentList(k)}.sqError-gMRAReg.Reg{i}{nodeList(k)}.sqError;
        end
    end
end

gMRAReg.Deltasq = gMRAReg.Deltasq/size(X,2);

gMRAReg.RegressionOpts = opts;

warning on;

gMRAReg.Timings.overall = cputime-gMRAReg.Timings.overall;
gMRAReg.Timings.overall_tic = toc;

return