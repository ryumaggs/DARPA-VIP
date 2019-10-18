function cX=Generate_Z_Manifold(NofPts, Dim, Opts)

%
% function cX=Generate_Z_Manifold(NofPts, Dim, Opts)
%
% Generate_Z_Manifold generates a Z shape manifold.
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
% Example: X = Generate_Z_Manifold(1000, 2, struct('DensityType', 'uniform'));
%          X = Generate_Z_Manifold(1000);
%

%% Setup parameters
if nargin < 2
    Dim = 2;
end

if nargin < 3
   Opts=[];
end

if ~isfield(Opts, 'DensityType')
    Opts.DensityType = 'uniform';
end
sigma = 0.5;

%% Generate data
%%%%%%%%% x1 x2 forms Z manifold x3,x4... uniform in [0,1)
cX     =zeros(NofPts,Dim+1);
NoPts1 = round(NofPts/(2+sqrt(2)));
NoPts3 = NoPts1;
NoPts2 = NofPts-2*NoPts1;
if strcmp(Opts.DensityType,'uniform') 
    dx1    = 2*(rand(NoPts1,1)-0.5);  
    dx1    = sort(dx1,'ascend');
    dx2    = 2*(rand(NoPts2,1)-0.5);  
    dx2    = sort(dx2,'descend');
    dx3    = 2*(rand(NoPts3,1)-0.5);  
    dx3    = sort(dx3,'ascend'); 
elseif strcmp(Opts.DensityType,'non-uniform') 
    dx1      = random(truncate(makedist('Normal','mu',1,'sigma',sigma),-1,1),NoPts1,1);
    dx1      = sort(dx1,'ascend');
    NoPts21  = round(NoPts2/2);
    dx21     = random(truncate(makedist('Normal','mu',1,'sigma',sigma),-1,1),NoPts21,1);
    NoPts22  = NoPts2-NoPts21;
    dx22     = random(truncate(makedist('Normal','mu',-1,'sigma',sigma),-1,1),NoPts22,1);
    dx2      = [dx21 ; dx22];
    dx2      = sort(dx2,'descend');
    dx3      = random(truncate(makedist('Normal','mu',-1,'sigma',sigma),-1,1),NoPts3,1);
    dx3      = sort(dx3,'ascend');
end
X3     = rand(NofPts, Dim-1);
cX(:,1:2)   = [dx1 , ones(NoPts1,1) ; dx2 , dx2 ; dx3 , -ones(NoPts3,1)] ;
cX(:,3:end) = X3;

    
return;

