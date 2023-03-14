% computes the euclidean projection onto the probability simplex
% https://arxiv.org/pdf/1309.1541.pdf
function X = SimplexProj(Y)

Y = Y';

[N,D] = size(Y);
X = sort(Y,2,'descend');
Xtmp = (cumsum(X,2)-1)*diag(sparse(1./(1:D)));
X = max(bsxfun(@minus,Y,Xtmp(sub2ind([N,D],(1:N),sum(X>Xtmp,2)))),0);

X = X';

end