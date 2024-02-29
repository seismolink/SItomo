function [L,fc] = likelyhood_smooth(data,indx,si,sig,A0,B0,C0,smoothpar)
% Copyright 2024 F.Link and M.D.Long 

M = length(si);

sig1 = A0.*cos(B0).*cos(C0);
sig2 = A0.*sin(B0).*cos(C0);
sig3 = A0.*sin(C0);
L1 = zeros(size(sig1));
L1(2:end-1,:,2:end-1) = sig1(1:end-2,:,2:end-1)+sig1(3:end,:,2:end-1)+sig1(2:end-1,:,1:end-2)+sig1(2:end-1,:,3:end)-4.*sig1(2:end-1,:,2:end-1);
L2 = zeros(size(sig2));
L2(2:end-1,:,2:end-1) = sig2(1:end-2,:,2:end-1)+sig2(3:end,:,2:end-1)+sig2(2:end-1,:,1:end-2)+sig2(2:end-1,:,3:end)-4.*sig2(2:end-1,:,2:end-1);
L3 = zeros(size(sig2));
L3(2:end-1,:,2:end-1) = sig3(1:end-2,:,2:end-1)+sig3(3:end,:,2:end-1)+sig3(2:end-1,:,1:end-2)+sig3(2:end-1,:,3:end)-4.*sig3(2:end-1,:,2:end-1);
mf2 = 1/2.*sum((L1(:)).^2+(L2(:)).^2+(L3(:)).^2)*length(si)/4;

fc = ((sum((data.si(indx)-si).^2))+mf2.*smoothpar.^2)./sig./sig./M;

L = exp(-fc/2);
end