function cX=Generate_CurvedZ_Manifold(NofPts, Dim, Opts)

%
% function cX=Generate_CurvedZ_Manifold(NofPts, Dim, Opts)
%
% Generate_CurvedZ_Manifold generates a curved Z shape manifold.
% 
% IN:
%    NofPts     : the number of points in the manifold generated
%    [Dim]      : the dimension of the manifold, if Dim=1, a curve, if Dim=2, a surface. default = 2
%    [Opts]     : structure containing the following fields:
%                   [DensityType] : 'uniform': uniform distribution 
%                                   'non-uniform': higher density at the
%                                   corners
%                                   default = non-uniform
%
% OUT:
%     cX: NofPtsxDim array, if PtsType is mesh, not exactly. 
%
% Example: X = Generate_CurvedZ_Manifold(1000, 2, struct('DensityType', 'non-uniform'));
%          X = Generate_CurvedZ_Manifold(1000);
%


% Setup parameters
if nargin < 2
    Dim = 2;
end

if nargin < 3
   Opts=[];
end

R = 5;  % curvature: 1/R


if ~isfield(Opts, 'DensityType')
    Opts.DensityType = 'uniform';
end
sigma = 0.25;



%%%%%%%%% x1 x2 forms Z manifold x3,x4... uniform in [0,1)
cX     =zeros(NofPts,Dim+1);
NoPts1 = round(NofPts/(2+sqrt(2)));
NoPts2 = NofPts-2*NoPts1;
% 
if strcmp(Opts.DensityType,'uniform') 
    alpha  = acos(1/R);
    dt     = (pi-2*alpha)*rand(NoPts1,1)+alpha;
    dt     = sort(dt,'descend');
    dx11   = R*cos(dt);
    dx12   = 1-sqrt(R^2-1)+R*sin(dt);
    %
    beta   = asin(sqrt(2)/R);
    dt     = 2*beta*rand(NoPts2,1)+(3*pi/4-beta);
    dt     = sort(dt,'ascend');
    dx21   = sqrt(R^2-2)*cos(-pi/4)+R*cos(dt);
    dx22   = sqrt(R^2-2)*sin(-pi/4)+R*sin(dt);
    %
    dt     = (pi-2*alpha)*rand(NoPts1,1)+alpha;
    dt     = sort(dt,'descend');
    dx31   = R*cos(dt);
    dx32   = 1-sqrt(R^2-1)+R*sin(dt);
elseif strcmp(Opts.DensityType,'non-uniform') 
    %
    alpha  = acos(1/R);
    dt     = (pi-2*alpha)*random(truncate(makedist('Normal','mu',0,'sigma',sigma),0,1),NoPts1,1)+alpha;
    dt     = sort(dt,'descend');
    dx11   = R*cos(dt);
    dx12   = 1-sqrt(R^2-1)+R*sin(dt);
    %
    NoPts21= round(NoPts2/2);
    NoPts22= NoPts2-NoPts21;
    beta   = asin(sqrt(2)/R);
    dt     = [2*beta*random(truncate(makedist('Normal','mu',1,'sigma',sigma),0,1),NoPts21,1)+(3*pi/4-beta);...
              2*beta*random(truncate(makedist('Normal','mu',0,'sigma',sigma),0,1),NoPts22,1)+(3*pi/4-beta)];
    dt     = sort(dt,'ascend');
    dx21   = sqrt(R^2-2)*cos(-pi/4)+R*cos(dt);
    dx22   = sqrt(R^2-2)*sin(-pi/4)+R*sin(dt);
    %
    dt     = (pi-2*alpha)*random(truncate(makedist('Normal','mu',1,'sigma',sigma),0,1),NoPts1,1)+alpha;
    dt     = sort(dt,'descend');
    dx31   = R*cos(dt);
    dx32   = 1-sqrt(R^2-1)+R*sin(dt);
    %
end

X3          = rand(NofPts, Dim-1);
cX(:,1:2)   = [dx11 , dx12 ; dx21 , dx22 ; dx31 , -dx32];
cX(:,3:end) = X3;
    
return;

