function [M, N, v] = spmtib(w, u)
% w has to be column vector

%w = w(end:-1:1);
%u = u(:, end:-1:1);

[p, n] = size(u);

[Mu, Nu] = sptrillrnk(u);
dw = spdiags(w, 0, n, n);
dcw = spdiags(sqrt(1-w.*conj(w)), 0, n, n);
Du = Nu - Mu;
M = (Nu + Du * conj(dw)) / dcw;
N = (Nu * dw + Du) / dcw;
v = zeros(n, p);
v(1:p, 1:p) = u(:,1:p)';
end