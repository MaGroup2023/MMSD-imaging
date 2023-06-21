function kk = k2kk(sz,k)
[i,j] = ind2sub(sz,k);
kk = sub2ind([sz(2),sz(1)],j,i);
end