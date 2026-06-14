function [M,N] = sptrillrnk(u)
[p, n] = size(u);
M = eye(n);
for i=1:n-p
    M(p+i, i:p+i-1) = -(u(:,i:i+p-1) \ u(:,i+p))';
end
M = sparse(M);
N = M * tril(u' * u);
N = sparse(N);
end