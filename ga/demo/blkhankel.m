function H = blkhankel(h)
% first block row to block Hankel
[di, do, l] = size(h);
k = floor((l + 1) / 2);
H = zeros(di * k, do * k);
for i = 1:k,
    for j = 1:k,
        H((i-1)*di+1:i*di, (j-1)*do+1:j*do) = h(:, :, i+j-1);
    end
end


% function T = blkhankel(v, m, n)
% 
%     % construct the (m,n)-block Hankel matrix with first n columns given by v
% 
%     p = length(v)/m;
% 
%     T = zeros(p*m, p*n);
%     T(:,1:n) = v;
%     sT2 = n;
%     j = [(m+1:m*p)'; (m*p+1)*ones(m,1)];
%     while 2*sT2 < n*p,
%         T(j<=m*p, sT2+(1:sT2)) = T(j(j<=m*p),1:sT2);
%         j(j<=m*p) = j(j(j<=m*p));
%         sT2 = sT2 * 2;
%     end
% 
%     T(j<=m*p, sT2+1:n*p) = T(j(j<=m*p), 1:(n*p-sT2));
% 
% end