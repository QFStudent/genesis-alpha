function y = hankelmul(h, x)
n = floor((length(h)+1)/2);
h = h(1:2*n-1);
hpad = [flipud(h); zeros(n-1,1)];
xpad = [x;zeros(2*(n-1), 1)];
ypad = ifft(fft(hpad) .* fft(xpad));
y = ypad(2*n-1:-1:n);
end