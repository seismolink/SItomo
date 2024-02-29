function Kt = transform_kernel3D(K,sets)
% Copyright 2024 F.Link and M.D.Long 

M = size(K);
Kt = zeros(size(K));
for i = 1:M(1)
    for j = 1:M(2)
        for k = 1:M(3)
            elem = K(sets(i,j,k).ix);
            Kt(i,j,k) = median(elem);
        end
    end
end
end