%% Demo of clustering
clear all;

pExampleNames   = { 
    'three gaussians','half moon','MNIST digits'};

fprintf('\n\n Select example to run:\n');
for k = 1:length(pExampleNames),
    fprintf('\n [%d] %s',k,pExampleNames{k});
end;
fprintf('\n\n  ');

if (~exist('pExampleIdx') || isempty(pExampleIdx) || pExampleIdx==0)
    pExampleIdx          = input('Pick an example to run:           ');
end

%% Set parameters for constructing the graph
EstDimOpts = struct('NumberOfTrials',15,'verbose',0,'MAXDIM',100,'MAXAMBDIM',100,'Ptwise',false,'NetsOpts',[],'UseSmoothedS',false, 'EnlargeScales',true );

Npts = 40000;

%% Generate data
switch pExampleIdx 
    case 1
        ncluster=3; r= 0.5;
        X1      = (mvnrnd([2,2],[r,0;0,r],500))';        Label1  = 1*ones(1,500);
        X2      = (mvnrnd([2,-2],[r,0;0,r],500))';       Label2  = 2*ones(1,500);
        X3      = (mvnrnd([-2*sqrt(3),0],[r,0;0,r],500))';       Label3  = 3*ones(1,500);
        X       = [X1 X2 X3];  Labels   = [Label1 Label2 Label3];
        ranper  = randperm(1500);
        X       = X(:,ranper); Labels = Labels(ranper);   
        clearvars X1 X2 X3 Label1 Label2 Label3 ranper
    case 2
        ncluster=2; noise=0.01;
        width   = 0.025*(rand(1,500)-0.5); center1 = [0 0]; center2=[1 0.5];
        theta1  = pi/2+pi*(rand(1,500)-0.5); % uniform %normrnd(pi/2,pi/7,1,500); % gaussian    
        X1(1,:) = center1(1)+(1+width).*cos(theta1); X1(2,:)=center1(2)+(1+width).*sin(theta1);   Label1 = 1*ones(1,500);
        theta2  = 3*pi/2+pi*(rand(1,500)-0.5);  %normrnd(3*pi/2,pi/7,1,500);  
        X2(1,:) = center2(1)+(1+width).*cos(theta2);  X2(2,:)=center2(2)+(1+width).*sin(theta2);  Label2 = 2*ones(1,500);
        X       = [X1 X2];  Labels   = [Label1 Label2];
        ranper  = randperm(1000);
        X       = X(:,ranper); Labels = Labels(ranper);  
        X       = X +  (mvnrnd([0,0],[noise,0;0,noise],1000))'; % add noise
        clearvars X1 X2 Label1 Label2 ranper center1 center2
    case 3
         % Loading data
         ncluster = 10;
         load TrainImages;  load TrainImageLabels;
         X=(reshape(TrainImages, [size(TrainImages, 1), size(TrainImages,2)*size(TrainImages,3), 1]))';
         X=double(X(:,1:end)); Labels=Labels(1:end);
end
        

% Set default options for constructing the graph
GraphDiffOpts = struct('Normalization','smarkov', 'Epsilon',1,'kNN', 45,'kNNAutotune', 10,'kEigenVecs', 35,'Symmetrization', 'W+Wt','DontReturnDistInfo', 1 );
Data.G        = GraphDiffusion(X, 0, GraphDiffOpts); 

% GraphDiffOpts = struct('Normalization','smarkov', 'Epsilon',sqrt(0.1),'kEigenVecs', 35,'Symmetrization', 'W+Wt','DontReturnDistInfo', 1 );
% Data.G        = GraphDiffusion(X, 0.2, GraphDiffOpts); 


% Original data colored by values of eigenfunctions
scrsz = get(groot,'ScreenSize');
figure('Position',[scrsz(3)*1/8,scrsz(4)*1/8,scrsz(3)*3/4,scrsz(4)*3/4]);
if pExampleIdx<3
    for k = 1:15
        subplot(3,5,k); scatter(X(1,:),X(2,:),20,Data.G.EigenVecs(:,k),'filled');    colorbar;
        title(sprintf('Eigenfunction %d on the data',k));
    end
elseif pExampleIdx == 3
   for k = 1:15
        subplot(3,5,k); imagesc(squeeze(TrainImages(k,:,:)));  
        title(sprintf('digits'));
   end
end

%% Cluster visulization

if ncluster == 2; markers = 'xo'; elseif ncluster == 3; markers = 'xo+';elseif ncluster == 10; markers = 'o*svhp<xh^>'; end
colors = hsv(ncluster);
figure('Position',[scrsz(3)*1/8,scrsz(4)*1/8,scrsz(3)*3/4,scrsz(4)*3/4]);
if pExampleIdx<3
    subplot(1,3,1); gscatter(X(1,:),X(2,:),Labels,colors,markers);
else
    [Q,~] = qr(randn(size(X,1)));
    QX    = Q*X;
    subplot(1,3,1); gscatter(QX(1,:),QX(2,:),Labels,colors,markers);
end
xlabel('x_1','Interpreter','tex');ylabel('x_2','Interpreter','tex');
title('Original data, first 2 coords');
subplot(1,3,2); gscatter(Data.G.EigenVecs(:,2),Data.G.EigenVecs(:,3),Labels,colors,markers);
xlabel('\phi_2','Interpreter','tex');ylabel('\phi_3','Interpreter','tex');
title('Diffusion embedding (2,3)');
subplot(1,3,3); gscatter(Data.G.EigenVecs(:,4),Data.G.EigenVecs(:,5),Labels,colors,markers);
xlabel('\phi_4','Interpreter','tex');ylabel('\phi_5','Interpreter','tex');
title('Diffusion embedding (4,5)');

fprintf('\n');

%% K-Means
if pExampleIdx==3
    kmeansLabels    = kmeans(Data.G.EigenVecs(:,1:11),10,'Replicates',10);
    [z,match]       = match_clusters(kmeansLabels,Labels);
    confmat = z;
    for k = 1:10
        confmat(:,k) = z(:,match(2,k));
    end
    fprintf('\n error rate=%.2f',1-sum(diag(confmat))/size(X,2));
    
end

        