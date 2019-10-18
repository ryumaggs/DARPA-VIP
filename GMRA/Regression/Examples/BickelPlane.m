clear all
close all
display('===============================')

%% Generate points from the plane x3=(x1+x2) wthin the cylinder
%% x1^2+x2^2<=25
D        = 50;
Nsample  = 1000;

X = zeros(D,Nsample);
for iter = 1 : Nsample
    x1   = (rand-0.5)*10;
    x2   = (rand-0.5)*10;
    while abs(x1)^2+abs(x2)^2>25
        x1   = (rand-0.5)*10;
        x2   = (rand-0.5)*10;
    end
    x3 = (x1 + x2);
    X(1:3,iter) = [x1 ; x2 ; x3];
end

figure
scatter3(X(1,:),X(2,:),X(3,:))


%% Regression function
ALs  = [ones(Nsample,1)  X.'];
%beta = rand(D+1,1).*sign(rand(D+1,1)-0.5);
load Beta.mat
Y    = ALs*beta;

display('---------------------------')
fprintf('exact coefficient = \n')
(beta(1:3)).'

%% Solve Pickel's LS --- to invert AtA !!!AtA only has rank 3
% AtA        = ALs.'*ALs/Nsample;
% [U EigAtA] = eig(AtA);
% % Fill zero with eps
% eps = 10^(-12);
% EigAtA = diag(EigAtA);
% EigAtA = max(EigAtA,eps);
% EigAtA = diag(EigAtA);
% 
% figure
% plot(sort(diag(EigAtA),'descend'),'ro-','LineWidth',2)
% axis tight
% title('eigenvalue in the noiseless case')
% 
% beta_rec = U*(EigAtA\(U'*(ALs.'*Y/Nsample)));
% 
% fX = ALs*beta_rec;       % regression at X
% ErrorPoint = abs(fX-Y);  % pointwise error
% 
display('---------------------------')
% display
fprintf('Bickel: Noise free case, rec coefficient = \n')
(beta_rec(1:3)).'

fprintf('relative L2 error =%6.4f percent\n',100*norm(ErrorPoint)/norm(Y))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Y is noisy
noiseY   = 50;  % noise percentage s.t. ||E||/||Y||=noise/100
sigma    = noiseY/100*norm(Y)/sqrt(Nsample);
Y_E      = normrnd(0,sigma,Nsample,1);
Y_noisy  = Y + Y_E;

%beta_rec = U*(EigAtA\(U'*(ALs.'*Y_noisy/Nsample)));

fX = ALs*beta_rec; % regression at Xi
ErrorPoint = abs(fX-Y);  % pointwise error

display('---------------------------')
% display
fprintf('Bickel: Y Noise = %6.0f percent, rec coefficient = \n',noiseY)
(beta_rec(1:3)).'

fprintf('relative L2 error =%6.4f percent\n',100*norm(ErrorPoint)/norm(Y))


% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display('---------------------------')
% 
% %% Add small gaussian noise on X
% noiseX     = 50;
% sigma      = 5*noiseX/100/sqrt(D);
% X_E        = normrnd(0,sigma,D,Nsample);
% X_noisy    = X + X_E; 
% ALs_noisy  = [ones(Nsample,1)  X_noisy.'];
% 
% AtA_noisy     = ALs_noisy.'*ALs_noisy/Nsample;
% [U_noisy EigAtA_noisy] = eig(AtA_noisy);
% 
% figure
% scatter3(X_noisy(1,:),X_noisy(2,:),X_noisy(3,:))
% 
% figure
% plot(sort(diag(EigAtA_noisy),'descend'),'ro-','LineWidth',2)
% axis tight
% title('eigenvalue in the presence of gaussian noise')
% 
% %% recovery with noise ! Noise is good for Bickel's method, matrix becomes invertible
% 
% beta_rec = U_noisy*(EigAtA_noisy\(U_noisy'*(ALs_noisy.'*Y/Nsample)));
% 
% fX = ALs_noisy*beta_rec; % regression at Xi
% ErrorPoint = abs(fX-Y);  % pointwise error
% 
% % display coefficent
% fprintf('Bickel: X noise = %6.0f\n',noiseX)
% fprintf('rec coefficient = \n')
% (beta_rec(1:3)).'
% 
% fprintf('relative L2 error =%6.4f percent\n',100*norm(ErrorPoint)/norm(Y))
% 
