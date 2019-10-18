% generate some data
N = 10000;
d = 2;
D = 3;
noise_level = 0;

XName           = 'D-Sphere'; XNickName = 'S';
XOpts           = struct('NumberOfPoints',N,'Dim',d,'EmbedDim',D,'NoiseType','Gaussian','NoiseParam',noise_level);
[X,~,OutParams] = GenerateDataSets( XName, XOpts);
OutParams.GraphDiffOpts.Display=0;

% split X into test and train sets
I           = randperm(N);
test_idx    = I(1:floor(.1*N));
train_idx   = I((1+floor(.1*N)):N);

% compute diffusion eigenvectors on training and full data (for comparison)
OutParams.GraphDiffOpts.kNN = 50;
Gtrain      = GraphDiffusion(X(:,train_idx),0,OutParams.GraphDiffOpts);
G           = GraphDiffusion(X,0,OutParams.GraphDiffOpts);

% run Nystrom extension 
Q = Nystrom(X(:,train_idx),X(:,test_idx),Gtrain,OutParams.GraphDiffOpts.kNN);

% compare with true diffusion coordinates
for m = 1:OutParams.GraphDiffOpts.kEigenVecs
    e(m) = min(norm(G.EigenVecs(test_idx,m)-Q(:,m)),norm(G.EigenVecs(test_idx,m)+Q(:,m))) / norm(G.EigenVecs(test_idx,m));
end

% plot results
figure;plot(e);
figure;imagesc(abs(normcols(G.EigenVecs(test_idx,:))'*normcols(Q)));colorbar;
