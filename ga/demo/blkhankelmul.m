function y = blkhankelmul(x, h)
% block Hankel multiplication
% H is first 'block row' of Hankel
[di, do, l] = size(h);
k = floor((l + 1) / 2);
y = 0 * x;
for i = 1:di,
    for j = 1:do,
        y(i:di:end) = y(i:di:end) + hankelmul(reshape(h(i, j, :), l, 1), x(j:do:end));
    end
end




