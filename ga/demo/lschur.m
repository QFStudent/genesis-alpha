function [L, Q] = lschur(A)
[U, R] = schur(A);
L = R(end:-1:1, end:-1:1);
Q = U(:, end:-1:1);
