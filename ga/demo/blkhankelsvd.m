function [svals, v] = blkhankelsvd(h, ord)
% problematic... for non-symmetric matrices, magnitude of eigenvalues is
% not the same as singular values
[di, do, l] = size(h);
k = ceil(l/2);
opts.disp = 0; opts.issym = false; opts.isreal = true; opts.tol = 10*eps;
[v, l] = eigs(@blkhankelmul, k*do, ord, 'LM', opts, h);
[~, j] = sort(abs(diag(l)), 'descend');
svals = abs(diag(l(j,j)));
[v,~] = qr(v(1:k*do, j(1:ord)), 0);
end