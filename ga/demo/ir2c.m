function c = ir2c(h, w, u)
[di, ~, l] = size(h);
n = length(w);
c = zeros(di, n);
[M, N, v] = spmtib(w, u);
K = M \ v;
for i=1:l,
    c = c + h(:, :, i) * K';
    K = M \ (N * K);
end