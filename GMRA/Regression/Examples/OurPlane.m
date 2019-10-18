close all
display('===============================')


%% SVD
Center    = sum(X,2)/Nsample;
Xcenter   = bsxfun(@minus,X,Center);
[U Xeig V] = svd(Xcenter);

figure
plot(diag(Xeig),'b*-','LineWidth',2)
title('singular value of X')

d     = 2;          % instrinsic dimension
U     = U(:,1:d);     % subspace
Xproj = U'*Xcenter;

ALs   = [ones(Nsample,1)  Xproj.'];


%% Noise free: solve our LS
AtA        = ALs.'*ALs/Nsample;
[U EigAtA] = eig(AtA);

figure
plot(sort(diag(EigAtA),'descend'),'ro-','LineWidth',2)
axis tight
title('eigenvalue in the noiseless case')


beta_rec = U*(EigAtA\(U'*(ALs.'*Y/Nsample)));

fX = ALs*beta_rec; % regression at Xi
ErrorPoint = abs(fX-Y);  % pointwise error

display('---------------------------')
% display
fprintf('SVD: Noise free case\n')

fprintf('relative L2 error =%6.4f percent\n',100*norm(ErrorPoint)/norm(Y))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Y is noisy
beta_rec = U*(EigAtA\(U'*(ALs.'*Y_noisy/Nsample)));

fX = ALs*beta_rec; % regression at Xi
ErrorPoint = abs(fX-Y);  % pointwise error

display('---------------------------')
% display
fprintf('SVD: Y Noise = %6.0f percent, rec coefficient = \n',noiseY)

fprintf('relative L2 error =%6.4f percent\n',100*norm(ErrorPoint)/norm(Y))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
display('---------------------------')

%% Add small gaussian noise on X -> X_noisy
% SVD
Center_noisy    = sum(X_noisy,2)/Nsample;
Xcenter_noisy   = bsxfun(@minus,X_noisy,Center_noisy);
[U_noisy Xeig_noisy V_noisy] = svd(Xcenter_noisy);

figure
plot(diag(Xeig_noisy),'b*-','LineWidth',2)
title('singular value of X_{noisy}')

d        = 2;          % instrinsic dimension
U_noisy  = U_noisy(:,1:d);     % subspace
Xproj_noisy = U_noisy'*Xcenter_noisy;

ALs_noisy   = [ones(Nsample,1)  Xproj_noisy.'];

AtA_noisy     = ALs_noisy.'*ALs_noisy/Nsample;
[U_noisy EigAtA_noisy] = eig(AtA_noisy);

figure
plot(sort(diag(EigAtA_noisy),'descend'),'ro-','LineWidth',2)
axis tight
title('eigenvalue in the presence of gaussian noise')

beta_rec = U_noisy*(EigAtA_noisy\(U_noisy'*(ALs_noisy.'*Y/Nsample)));

fX = ALs_noisy*beta_rec; % regression at X
ErrorPoint = abs(fX-Y);  % pointwise error

% display coefficent
fprintf('SVD: X noise = %6.0f\n',noiseX)

fprintf('relative L2 error =%6.4f percent\n',100*norm(ErrorPoint)/norm(Y))

