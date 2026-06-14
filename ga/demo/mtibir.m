function ir = mtibir(lambda, u, C, m)
% this is not the optimized version
[M, N, v] = spmtib(lambda, u);
[di, d] = size(u);
do = size(C, 1);
J = zeros(d, di, m);
J(:,:,1) = M \ v;
for i = 2:m
    J(:,:,i) = M \ (N * J(:,:,i-1));
end
ir = zeros(do, di, m);
for i = 1:m
    ir(:,:,i) = C * J(:,:,i);
end