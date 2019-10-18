% Generate samples
fprintf('\n Generate Samples...');
[X, GMRAopts, imgOpts, X0, Labels, OutParams] = GenerateData_and_SetParameters(Params.ExampleName,Params.N,Params.d);

if isempty(X0), X0 = X; end;

X       = single(X);        
X0      = single(X0);

% Add Gaussian noise
if Params.Sigma > 0
    X       = X+Params.Sigma*randn(size(X))/sqrt(size(X,1));
end

% Create train and test data
TrainAndTestOpts.TestSize = floor(Params.TestPercent*size(X,2)/100);
if TrainAndTestOpts.TestSize > 0
    TrainAndTestOpts.TrainIdxs                  = randperm(size(X,2),size(X,2)-TrainAndTestOpts.TestSize);
    TrainAndTestOpts.TestIdxs                   = setdiff(1:size(X,2),TrainAndTestOpts.TrainIdxs);
    TrainAndTestOpts.TrainLabels(size(X,2),1)   = false;
    TrainAndTestOpts.TrainLabels(TrainAndTestOpts.TrainIdxs) = true;
    X_train                                     = X(:,TrainAndTestOpts.TrainIdxs);
    X_test                                      = X(:,TrainAndTestOpts.TestIdxs);
    X0_train                                    = X0(:,TrainAndTestOpts.TrainIdxs);
    X0_test                                     = X0(:,TrainAndTestOpts.TestIdxs);
else
    TrainAndTestOpts.TrainIdxs                  = 1:size(X,2);
    TrainAndTestOpts.TestIdxs                   = [];
    TrainAndTestOpts.TrainLabels(1:size(X,2),1) = true;
    X_train                                     = X;
    X_test                                      = [];
end
