function [w, u, C] = msvdreduce(h, order)
[~, do, ~] = size(h);
[~, ~, V] = svd(blkhankel(h));
V = V(:, 1:order);
A = V(do+1:end, :)' * V(1:end-do, :);
[A, Q] = lschur(A);
B = Q' * V(1:do, :)';

u = zeros(do, order);
w = diag(A);
for i=1:order,
    u(:, i) = B(1, :)' / norm(B(1, :));
    B = B(2:end, :);
    B = B - (1 + 1/conj(w(i))) * (B * u(:, i)) * u(:, i)';
end
C = ir2c(h, w, u);
end

