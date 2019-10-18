function Y_test = NNclassify( X_train, Y_train, X_test, K, coverTree, isClassificationProblem )

if nargin<6, isClassificationProblem = true; end

% Find nearest neighbors of each test point
[~, kNNIdx]    = nrsearch(X_train, X_test, max(K), 0, struct('ReturnAsArrays',1,'NNInfo',struct('CoverTree',coverTree)));

% Find most frequent label among the nn's
lUniqueTrainingLabels = unique(Y_train);
Y_test = zeros(size(X_test,2),length(K));

if isClassificationProblem
    for l = 1:length(K)
        for k = 1:size(X_test,2)
            histk           = hist(Y_train(kNNIdx(k,1:K(l))),lUniqueTrainingLabels);
            [~,Y_test(k,l)] = max(histk);
        end
    end
else
    for l = 1:length(K)
        for k = 1:size(X_test,2)
            Y_test(k,l) = mean(Y_train(kNNIdx(k,1:K(l))));
        end
    end
end


if isClassificationProblem
    Y_test = lUniqueTrainingLabels(Y_test);
end

return