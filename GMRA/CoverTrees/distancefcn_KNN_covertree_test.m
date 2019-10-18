%% Test for distancefcn_KNN_covertree
% (c) Mauro Maggioni

INCLUDE_COVERTREE_DISTANCES

%% Euclidean distances
kNN = 10;

X_1 = randn(30,100); 
X_2 = randn(30,200); 
fprintf('\n Covertree nn''s...');
tic;
[D,idxs,dists,NNInfo] = distancefcn_KNN_covertree( X_1, X_2, [], COVERTREE_DISTANCE_EUCLIDEAN_ABS_VALUES, ...
                                struct('data_classname',COVERTREE_CLASSNAME_VECTOR,'ReturnD', false,'kNN',kNN,'NTHREADS',0) );
timings.distancefcn_KNN_covertree = toc;                            
fprintf('done.');

% Check correctness by comparing with brute-force computation
fprintf('\n All pairwise...');
tic
dists_bf                = pdist2(X_1',X_2');
timings.dists_bf        = toc;
[sorted,sorted_idxs]    = sort(dists_bf,1,'ascend');
sorted                  = sorted(1:kNN,:);
sorted_idxs             = int32(sorted_idxs(1:kNN,:));
fprintf('done.');

fprintf('\n Maximum difference (covertree vs. brute force) in index of nearest neighbors (should be 0):%d',max(max(sorted_idxs-idxs)));
fprintf('\n Maximum difference (covertree vs. brute force) in distances to nearest neighbors (should be O(eps)):%e',max(max(sorted-dists)));

% Go multi-threaded
fprintf('\n Covertree nn''s, multi-threaded...');
tic
[D,idxs,dists,NNInfo] = distancefcn_KNN_covertree( X_1, X_2, [], COVERTREE_DISTANCE_EUCLIDEAN, ...
                                struct('data_classname',COVERTREE_CLASSNAME_VECTOR,'ReturnD', false,'kNN',kNN,'NTHREADS',4) );
timings.distancefcn_KNN_covertree_mt = toc;
fprintf('done.');

fprintf('\n Maximum difference (multi-threaded covertree vs. brute force) in index of nearest neighbors (should be 0):%d',max(max(sorted_idxs-idxs)));
fprintf('\n Maximum difference (multi-threaded  covertree vs. brute force) in distances to nearest neighbors (should be O(eps)):%e',max(max(sorted-dists)));

timings

%% Molecular distances - just to show to call the function when using a different metric
X_1 = randn(30,100); 
X_2 = randn(30,200); 
[D,idxs,dists,NNInfo] = distancefcn_KNN_covertree( X_1, X_2, [], COVERTREE_DISTANCE_RMSD, ...
                                struct('data_classname',COVERTREE_CLASSNAME_MOLECULARSTATE,'ReturnD', true,'kNN',kNN) );
                            
fprintf('\n');                            
