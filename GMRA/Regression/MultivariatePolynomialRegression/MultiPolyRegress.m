function reg = MultiPolyRegress(X,Y,PW,varargin)

%
%   reg = MultiPolyRegress(Data,R,PW) performs multi-variable polynomial
%   regression analysis on row stacked dimensional data matrix Data. Data is
%   an m-by-n (m>n) matrix where n is the number of data points and m is the number of
%   independent variables. R is the m-by-1 response column vector and PW is the degree
%   of the polynomial fit.
%
%   reg = MultiPolyRegress(...,PV) restricts individual dimensions of
%   Data to particular powers PV in the polynomial expansion. PV is an
%   m-by-1 vector. A PV of [2 1] would limit a 2-dimensional 2nd degree polynomial to
%	the terms that have x^2, x and y, eliminating the terms with y^2.
%
%
%   reg = MultiPolyRegress(...,'range') adjusts the normalization of
%   goodness of fit measures: mean of absolute error (mae) and standard deviation
%   of absolute error. i.e. by default, mae is defined mean(abs(y-yhat)./y),
%   however, when this switch is used, the definition is changed to
%   mean(abs(y-yhat)./range(y)). It is useful when your y vector (R in the
%   syntax of this code) contains values that are equal to or very close to
%   0.
%
%   reg is a struct with the following fields:
%          FitParameters: Section Header
%            PowerMatrix: A matrix that describes the powers for each term
%                         of the polynomial fit. It is useful for
%                         evaluating any future points with the calculated
%                         fit. Refer to the "Compose" section in the code
%                         on how to use it.
%                 Scores: Is a diagnostic reference. Displays the raw value
%                         of individual polynomial terms for each data
%                         point, before multiplication with coefficients.
%                         In other words, it is the matrix X you would have
%                         input in to the Statistical Toolbox function
%                         "regress".
%   PolynomialExpression: The expression for the fitted polynomial.
%           Coefficients: For the calculated fit.
%                   yhat: Estimated values by the fit.
%              Residuals: y-yhat or R-yhat in syntax of this code,
%          GoodnessOfFit: Section Header
%                RSquare: 1-SSE/TSE
%                    MAE: Normalized Mean of Absolute Error
%                 MAESTD: Standard Deviation of Absolute Error
%     LOOCVGoodnessOfFit: '-----------------'
%              CVRSquare: RSquare of LOOCV
%                  CVMAE: MAE of LOOCV
%               CVMAESTD: MAESTD of LOOCV
%
%   Copyright (c) 2015, Ahmet Cecen  -  All rights reserved.
%   Modified by Mauro Maggioni [in a non-backward-compatible way]
%

warning off 

% check if empty
if isempty(X), reg = struct('Coefficients',[],'sqError',0); return; end;

% Arrange Input Arguments
PV = repmat(PW,[1,size(X,1)]);
Reg = [];
NormalizationSwitch='1-to-1 (Default)';
if nargin > 3
    for ii=1:length(varargin)
        if strcmp(varargin{ii},'range'),        NormalizationSwitch='Range';    end
        if strcmp(varargin{ii},'reg'),          Reg = varargin{ii+1};           end
        if strcmp(varargin{ii},'PV'),           PV=varargin{ii+1};              end
    end
end

% Function Parameters
NData = size(X,2);
NVars = size(X,1);
RowMultiC = '1';
Lim = max(PV);

if isempty(Reg)
    if Lim>1,
        A = zeros(Lim^NVars,NVars);                                                                                                     % Initialize
        for ii=1:NVars                                                                                                                  % Create Colums Corresponding to Mathematical Base
            A(:,ii)=mod(floor((1:Lim^NVars)/Lim^(ii-1)),Lim);
        end
        
        A=fliplr(A); A=A(sum(A,2)<=Lim,:); Ab=diag(repmat(Lim,[1,NVars])); A=[A;Ab];                                                    % Flip - Reduce - Augment
        
        for ii=1:NVars                                                                                                                  % Degree Conditionals
            A=A(A(:,ii)<=PV(ii),:);
        end
        
        NLegend = size(A,1);                                                                                                            % Build Framework
        Scores = zeros(NData,NLegend);                                                                                                  % Allocate
        for ii=1:NVars
            RowMultiC=strcat(RowMultiC,['.*C(:,',num2str(ii),')']);
        end
        for ii=1:NData                                                                                                                  % Compose
            current=repmat(X(:,ii)',[NLegend,1]);
            C=current.^A; %#ok<NASGU>
            Scores(ii,:) = eval(RowMultiC);
        end
    elseif Lim==1
        A      = [0,0;1,0;0,1];
        Scores = [ones(1,NData);X]';
    elseif Lim==0
        A      = [0,0];
        Scores = ones(NData,1);
    end
else
    A = Reg.PowerMatrix;
    Scores = Reg.Scores;
end

[QQ,RR,perm] = qr(Scores,0);                                                                                                    % Use  QR to avoid explicit inversion and check rank.

p = sum(abs(diag(RR)) > size(Scores,2)*eps(RR(1)));

if true, %p < size(Scores,2)
    %warning('MultiPolyRegress:Rank Deficiency within Polynomial Terms!');
    RR = RR(1:p,1:p);
    QQ = QQ(:,1:p);
    perm = perm(1:p);
end

b       = zeros(size(Scores,2),size(Y,2));                                                                                              % Ordinary Least Squares Regression
if ~isempty(Y)
    b(perm,:)   = RR \ (QQ'*Y);
    yhat        = Scores*b;
    r           = Y-yhat;
    sqError     = sum(r.^2);
else
    yhat        = [];
    r           = [];
    sqError     = [];
end

% % Goodness of Fit
%r2 = 1 - (norm(r)).^2/norm(Y-mean(Y))^2;
% if strcmp(NormalizationSwitch,'Range')==1
%     mae = mean(abs(r./abs(max(R)-min(R))));
%     maestd = std(abs(r./abs(max(R)-min(R))));
% else
%     mae = mean(abs(r./R));
%     maestd = std(abs(r./R));
% end
%
% % Leave One Out Cross Validation
% H=QQ*QQ';
% rCV=r./(1-diag(H));
%
% % LOOCV Goodness of Fit
% CVr2 = 1 - (norm(rCV)).^2/norm(R-mean(R))^2;
% if strcmp(NormalizationSwitch,'Range')==1
%     CVmae = mean(abs(rCV./abs(max(R)-min(R))));
%     CVmaestd = std(abs(rCV./abs(max(R)-min(R))));
% else
%     CVmae = mean(abs(rCV./R));
%     CVmaestd = std(abs(rCV./R));
% end

% Construct Output
% reg = struct('FitParameters','-----------------','PowerMatrix',A,'Scores',Scores, ...
%     'Coefficients',b, 'yhat', yhat, 'Residuals', r, ...
%     'GoodnessOfFit','-----------------', 'RSquare', r2, 'MAE', mae, 'MAESTD', maestd, ...
%     'Normalization',NormalizationSwitch,'LOOCVGoodnessOfFit','-----------------', 'CVRSquare', ...
%     CVr2, 'CVMAE', CVmae, 'CVMAESTD', CVmaestd,'CVNormalization',NormalizationSwitch);

reg = struct('PowerMatrix',A,'Lim',Lim,'Scores',Scores, 'Coefficients',b, 'yhat', yhat, 'Residuals', r, 'sqError', sqError, ...
    'Normalization',NormalizationSwitch);%,'r2',r2);


return

