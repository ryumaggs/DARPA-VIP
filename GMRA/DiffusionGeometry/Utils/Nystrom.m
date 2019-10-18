function Y = Nystrom(X, Xq, G, kNN)
%
%
% function Y = Nystrom(X, Xq, G, kNN)
%
% NYSTROM Nystrom extension of diffusion coordinates to new data points
% input:     
%            X: original data points (D x N)
%            Xq: new data points (D x L)
%            G: graph structure returned by GraphDiffusion on X
%            kNN: number of nearest neighbors used in Nystrom extension
% output:
%            Y: approximate diffusion coordinates for Xq
%
% Author: David Lawlor (djl@math.duke.edu)
%   Mauro Maggioni (mauro.maggioni@jhu.edu)
%

[D,N] = size(X);
[D1,L] = size(Xq);

if (D~=D1)
    fprintf(1,'Error; X and Xq must have same number of rows\n');
    return
end

% calculate distances of new points to old
[count, idxs, dists] = nrsearch(X, Xq, kNN, 0, []);

% now construct the L x N weight matrix W (copied from GraphDiffusion)
rowidxs = zeros(sum(count), 1,'double');
colidxs = zeros(sum(count), 1,'double');
distsmat = zeros(sum(count), 1);
epsmat = zeros(length(count));
index = 1;
location = 1;
for j=1:length(count)
    colidxs(location:(location+count(j)-1))     = idxs{index};
    rowidxs(location:(location+count(j)-1))     = index;
    if isfield(G,'Autotune')
      epsmat(index) = mean(G.Autotune(idxs{index}));
    else
      epsmat(index) = G.Epsilon;
    end
    distsmat(location:(location+count(j)-1))    = dists{index}/epsmat(index);
    location                                    = location+count(j);
    index                                       = index+1;
end

W = sparse(rowidxs, colidxs,  exp(-distsmat.^2), L, N, sum(count));
D = sum(W,2);
W = bsxfun(@rdivide,W,D);
Y = W * G.EigenVecs;% * spdiags(1./G.EigenVals,0,length(G.EigenVals),length(G.EigenVals));

%DInvSqrt = zeros(size(D));
%lNonZeroIdxs = find(D~=0);
%DInvSqrt(lNonZeroIdxs) = D(lNonZeroIdxs).^(-1);
%DInvSqrt = sqrt(DInvSqrt);
%DInvSqrt = sparse(1:N, 1:N, DInvSqrt, N, N, N);
%T = DInvSqrt*W*DInvSqrt;
%T = (T'+T)/2;
%Y = T * G.EigenVecs * spdiags(1./G.EigenVals,0,length(G.EigenVals),length(G.EigenVals));

return
