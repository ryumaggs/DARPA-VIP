% Function to generate data from PumaDyn simulator.
% Based on (Ghahramani,1996).
% Requires Robotics Toolbox (Corke,1996).
%
% INPUT:  j     = number of free joints (from 1 to 6, default=3)
%         N     = number of samples (default=10000)
%         beta  = nonlinearity parameter (default=1)
%         gamma = unpredictability parameter (default=0)
%
% OUTPUT: X     = N-by-d matrix of input points (d=3j-1)
%         Y     = N=by-1 vector of output values
%
% REQUIREMENTS
%   puma560 from '96 version of http://petercorke.com/Robotics_Toolbox.html
%
% PumaDyn-8 dataset parameters: j = 3 (d=8), N = 8192
%                             _________________
%                            |  beta  |  gamma |
%                        nm  |  1.2   |  0.12  |
%                        nh  |  1.2   |  0.40  |
%                        fm  |  0.6   |  0.12  |
%                        fh  |  0.6   |  0.40  |
%                             _________________


function [X,Y] = puma(j,N,beta,gamma)

puma560; % read parameters of the Puma arm

if nargin == 0
    j = 3;
end

if nargin <=1
    N     = 10000;
end

if nargin <= 2
    beta = 1;
end

if nargin <= 3
    gamma = 0;
end

d = 3*j-1;

X = zeros(N,d);
Y = zeros(N,1);

for i=1:N;
    thetan  = pi*[(rand(j,1)-0.5)*beta;zeros(6-j,1)];      % generate positions
    thetadn = pi*[(rand(j,1)-0.5)*beta;zeros(6-j,1)];      % generate velocities
    torquen = [(rand(j-1,1)-0.5)*beta;zeros(6-j+1,1)];     % generate torques
    
    theta   = thetan  + randn(6,1)*gamma;                  % add Gaussian noise
    thetad  = thetadn + randn(6,1)*gamma;                  % to inputs
    torque  = torquen + randn(6,1)*gamma;
    
    thetadd = accel(p560,theta,thetad,torque);             % compute forward dynamics
    
    X(i,:) = [thetan(1:j)' thetadn(1:j)' torquen(1:j-1)']; % assemble input vector
    Y(i,:) = thetadd(3) + gamma*randn;                     % add noise to target
end
        
end