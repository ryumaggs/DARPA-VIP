% EXAMPLE of usage of covertree_expand
X = randn(2,10000);
opts = struct('theta',0.5,'numlevels',int32(1000),'minlevel',int32(0),'NTHREADS',int32(0),'BLOCKSIZE',int32(1024));
CoverTree = covertree( opts, X);
CoverTree_exp1 = covertree_expand( CoverTree, struct('ExtRange','max','MaxLevel',Inf) );
CoverTree_exp2 = covertree_expand( CoverTree, struct('ExtRange','min','MaxLevel',4) );
figure;treeplot(double(CoverTree.levels(:,2)+1)');
%figure;treeplot(double(CoverTree_exp1.levels(:,2)+1)');
figure;treeplot(double(CoverTree_exp2.levels(:,2)+1)');


