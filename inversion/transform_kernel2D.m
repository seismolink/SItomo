function Kt = transform_kernel2D(K,sets)
% Copyright 2024 F.Link and M.D.Long 

M = size(K);
K2 = squeeze(K);
for i = 1:M(1)
    for j = 1:M(3)
        elem = K2(sets(i,j).ix);
        Kn(i,j) = median(elem);
    end
end
Kt = reshape(Kn,M(1),M(2),M(3));
end