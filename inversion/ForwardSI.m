 function [si,vi] = ForwardSI(per,par,AZIMUTH,x_1,y_1,z_1,X,Y,Z,omega,alpha_0,beta_0,strength_in,phi_in,theta_in,vflag)
% Copyright 2024 F.Link and M.D.Long 

[M1,M2,M3] = size(strength_in);
temp = unique(X);
dx = temp(1)-temp(2);
temp = unique(Y);
dy = temp(1)-temp(2);
temp = unique(Z);
dz = temp(1)-temp(2);


% initialize modelspace
iota=sqrt(-1);

% Initial waveform
m=((omega.^2.*per.^2)./(4.*pi)).*exp(-((omega.^2.*per.^2)./(8.*pi.^2)));
dom = omega(2)-omega(1);
normf = sum(m.*dom);
m = m./normf;

%%%%%%%%%%%%%%%%%%%%
% START CALCULATIONS
%%%%%%%%%%%%%%%%%%%%

if vflag 
    inc = zeros(size(beta_0));
else
    inc = asind(par./1000.*beta_0);
end

% position vector
R=sqrt((X-x_1).^2+(Y-y_1).^2+(Z-z_1).^2);

% wavevector 
wv(1,:) = -sind(inc(:)).*sind(AZIMUTH+180);
wv(2,:) = sind(inc(:)).*cosd(AZIMUTH+180); % here a negative sign?
wv(3,:) = cosd(inc(:));

% Z = -Z;

% detour of scattered wave relative to direct wave
D = abs(reshape(dot(wv,[X(:)-x_1 Y(:)-y_1 Z(:)-z_1]',1),M1,M2,M3));

% Azimuth of wave and polarization
Alpha=AZIMUTH.*pi/180;
Alpha2 = AZIMUTH.*pi/180;

% computing the components of the unit vector
p1=(X-x_1)./R; 
p2=(Y-y_1)./R;
p3=(Z-z_1)./R;

%calculating the direction cosines
s1=-cos(-theta_in).*sin(phi_in);
s2=cos(-theta_in).*cos(phi_in);
s3=sin(-theta_in);

%Polarization vector
gg(1,:) = -cosd(inc(:)).*sind(AZIMUTH+180);
gg(2,:) = cosd(inc(:)).*cosd(AZIMUTH+180);
gg(3,:) = -sind(inc(:));

% phase, far field pre-factor, second derivative
M = iota.*exp((iota.*(R-D)./beta_0).*reshape(omega,1,1,1,[])).*(ones(M1,M2,M3).*reshape(m,1,1,1,[])).*dom;
sumsum = sum(omega.^2.*m.*dom);
fFF = reshape(omega.^3,1,1,1,[])./4./pi./beta_0.^2./R;
fMFF = sum(real(fFF.*M),4);
gFF=real(fMFF)./(sumsum);

fMF = -iota.*reshape(omega.^2,1,1,1,[])./4./pi./beta_0./R.^2;
fMMF = sum(real(fMF.*M),4);
gMF=real(fMMF)./(sumsum);

fNF = -1./4./pi./R.^3;
fMNF = fNF.*sum(real(M).*reshape(omega,1,1,1,[]),4);
gNF=real(fMNF)./(sumsum);

%%%% new until here

% Radiation pattern
wv1 = reshape(wv(1,:),M1,M2,M3);
wv2 = reshape(wv(2,:),M1,M2,M3);
wv3 = reshape(wv(3,:),M1,M2,M3);
gg1 = reshape(gg(1,:),M1,M2,M3);
gg2 = reshape(gg(2,:),M1,M2,M3);
gg3 = reshape(gg(3,:),M1,M2,M3);

% Near field
exR = -6.*(5.*p1.^3.*(s1.^2-1)+10.*p1.^2.*s1.*(p2.*s2+p3.*s3)+p1.*(5.*p2.^2.*(s2.^2-1)+10.*p3.*p2.*s2.*s3+5.*p3.^2.*(s3.^2-1)-3.*s1.^2-s2.^2-s3.^2+5)-2.*s1.*(p2.*s2+p3.*s3)).*(-gg1.*s1.*(s2.*wv2+s3.*wv3)+gg3.*((s1.^2-s3.^2).*wv3-s3.*(s1.*wv1+s2.*wv2))+gg2.*(s1.^2.*wv2-s1.*s2.*wv1-s2.*(s2.*wv2+s3.*wv3)));
eyR = -6.*(5.*p2.^3.*(s2.^2-1)+10.*p2.^2.*s2.*(p1.*s1+p3.*s3)+p2.*(5.*p1.^2.*(s1.^2-1)+10.*p3.*p1.*s1.*s3+5.*p3.^2.*(s3.^2-1)-s1.^2-3.*s2.^2-s3.^2+5)-2.*s2.*(p1.*s1+p3.*s3)).*(-gg1.*s1.*(s2.*wv2+s3.*wv3)+gg3.*((s1.^2-s3.^2).*wv3-s3.*(s1.*wv1+s2.*wv2))+gg2.*(s1.^2.*wv2-s1.*s2.*wv1-s2.*(s2.*wv2+s3.*wv3)));
gxR = 6.*(-2.*p1.*(5.*p3.^2-1).*(-gg1.*s1.*s2.*wv2+gg2.*((s1.^2-s2.^2).*wv2-s1.*s2.*wv1)+gg3.*(s1.^2+s3.^2).*wv3)+(5.*p1.^2-1).*p3.*(4.*gg2.*s1.*s3.*wv2-2.*s2.*s3.*(gg2.*wv1+gg1.*wv2)-2.*(s1.^2+s3.^2).*(gg3.*wv1+gg1.*wv3)-2.*s1.*s2.*(gg3.*wv2+gg2.*wv3))+(5.*p1.^2-1).*p2.*(-2.*(s1.^2+s2.^2).*(gg2.*wv1+gg1.*wv2)+4.*gg3.*s1.*s2.*wv3-2.*s2.*s3.*(gg3.*wv1+gg1.*wv3)-2.*s1.*s3.*(gg3.*wv2+gg2.*wv3))+5.*p1.*p2.*p3.*(4.*gg1.*s2.*s3.*wv1-2.*s1.*s3.*(gg2.*wv1+gg1.*wv2)-2.*s1.*s2.*(gg3.*wv1+gg1.*wv3)-2.*(s2.^2+s3.^2).*(gg3.*wv2+gg2.*wv3))-2.*p1.*(5.*p2.^2-1).*(gg2.*(s1.^2+s2.^2).*wv2-gg1.*s1.*s3.*wv3+gg3.*((s1.^2-s3.^2).*wv3-s1.*s3.*wv1))+2.*p1.*(5.*p1.^2-3).*(gg2.*((s1.^2+s2.^2).*wv2+s2.*s3.*wv3)+gg3.*(s2.*s3.*wv2+(s1.^2+s3.^2).*wv3)));
gyR = 6.*(-2.*p2.*(5.*p3.^2-1).*(-gg1.*s1.*s2.*wv2+gg2.*((s1.^2-s2.^2).*wv2-s1.*s2.*wv1)+gg3.*(s1.^2+s3.^2).*wv3)+2.*p2.*(5.*p2.^2-3).*(gg1.*(s1.^2+s2.^2).*wv1+gg3.*(s2.^2+s3.^2).*wv3+s1.*s3.*(gg3.*wv1+gg1.*wv3))+5.*p1.*p2.*p3.*(4.*gg2.*s1.*s3.*wv2-2.*s2.*s3.*(gg2.*wv1+gg1.*wv2)-2.*(s1.^2+s3.^2).*(gg3.*wv1+gg1.*wv3)-2.*s1.*s2.*(gg3.*wv2+gg2.*wv3))+p1.*(5.*p2.^2-1).*(-2.*(s1.^2+s2.^2).*(gg2.*wv1+gg1.*wv2)+4.*gg3.*s1.*s2.*wv3-2.*s2.*s3.*(gg3.*wv1+gg1.*wv3)-2.*s1.*s3.*(gg3.*wv2+gg2.*wv3))+(5.*p2.^2-1).*p3.*(4.*gg1.*s2.*s3.*wv1-2.*s1.*s3.*(gg2.*wv1+gg1.*wv2)-2.*s1.*s2.*(gg3.*wv1+gg1.*wv3)-2.*(s2.^2+s3.^2).*(gg3.*wv2+gg2.*wv3))+2.*(5.*p1.^2-1).*p2.*(gg2.*((s1.^2+s2.^2).*wv2+s2.*s3.*wv3)+gg3.*(s2.*s3.*wv2+(s1.^2+s3.^2).*wv3)));
eNFTxy = (cosd(inc)).*cos(Alpha2).*2.*alpha_0.^2./beta_0.^2.*exR+(cosd(inc)).*sin(Alpha).*2.*alpha_0.^2./beta_0.^2.*eyR;
gNFTxy = 2.*((cosd(inc)).*cos(Alpha2).*gxR+(cosd(inc)).*sin(Alpha2).*gyR);
% Mid field
exR = -2.*(6.*p1.^3.*(s1.^2-1)+12.*p1.^2.*s1.*(p2.*s2+p3.*s3)+p1.*(6.*p2.^2.*(s2.^2-1)+12.*p3.*p2.*s2.*s3+6.*p3.^2.*(s3.^2-1)-4.*s1.^2-s2.^2-s3.^2+6)-3.*s1.*(p2.*s2+p3.*s3)).*(-gg1.*s1.*(s2.*wv2+s3.*wv3)+gg3.*((s1.^2-s3.^2).*wv3-s3.*(s1.*wv1+s2.*wv2))+gg2.*(s1.^2.*wv2-s1.*s2.*wv1-s2.*(s2.*wv2+s3.*wv3)));
eyR = -2.*(6.*p2.^3.*(s2.^2-1)+12.*p2.^2.*s2.*(p1.*s1+p3.*s3)+p2.*(6.*p1.^2.*(s1.^2-1)+12.*p3.*p1.*s1.*s3+6.*p3.^2.*(s3.^2-1)-s1.^2-4.*s2.^2-s3.^2+6)-3.*s2.*(p1.*s1+p3.*s3)).*(-gg1.*s1.*(s2.*wv2+s3.*wv3)+gg3.*((s1.^2-s3.^2).*wv3-s3.*(s1.*wv1+s2.*wv2))+gg2.*(s1.^2.*wv2-s1.*s2.*wv1-s2.*(s2.*wv2+s3.*wv3)));
gxR = -4.*p1.*(6.*p3.^2-1).*(-gg1.*s1.*s2.*wv2+gg2.*((s1.^2-s2.^2).*wv2-s1.*s2.*wv1)+gg3.*(s1.^2+s3.^2).*wv3)+2.*(3.*p1.^2-1).*p3.*(4.*gg2.*s1.*s3.*wv2-2.*s2.*s3.*(gg2.*wv1+gg1.*wv2)-2.*(s1.^2+s3.^2).*(gg3.*wv1+gg1.*wv3)-2.*s1.*s2.*(gg3.*wv2+gg2.*wv3))+(6.*p1.^2-1).*p3.*(4.*gg2.*s1.*s3.*wv2-2.*s2.*s3.*(gg2.*wv1+gg1.*wv2)-2.*(s1.^2+s3.^2).*(gg3.*wv1+gg1.*wv3)-2.*s1.*s2.*(gg3.*wv2+gg2.*wv3))+2.*(3.*p1.^2-1).*p2.*(-2.*(s1.^2+s2.^2).*(gg2.*wv1+gg1.*wv2)+4.*gg3.*s1.*s2.*wv3-2.*s2.*s3.*(gg3.*wv1+gg1.*wv3)-2.*s1.*s3.*(gg3.*wv2+gg2.*wv3))+(6.*p1.^2-1).*p2.*(-2.*(s1.^2+s2.^2).*(gg2.*wv1+gg1.*wv2)+4.*gg3.*s1.*s2.*wv3-2.*s2.*s3.*(gg3.*wv1+gg1.*wv3)-2.*s1.*s3.*(gg3.*wv2+gg2.*wv3))+12.*p1.*p2.*p3.*(4.*gg1.*s2.*s3.*wv1-2.*s1.*s3.*(gg2.*wv1+gg1.*wv2)-2.*s1.*s2.*(gg3.*wv1+gg1.*wv3)-2.*(s2.^2+s3.^2).*(gg3.*wv2+gg2.*wv3))-4.*p1.*(6.*p2.^2-1).*(gg2.*(s1.^2+s2.^2).*wv2-gg1.*s1.*s3.*wv3+gg3.*((s1.^2-s3.^2).*wv3-s1.*s3.*wv1))+8.*p1.*(3.*p1.^2-2).*(gg2.*((s1.^2+s2.^2).*wv2+s2.*s3.*wv3)+gg3.*(s2.*s3.*wv2+(s1.^2+s3.^2).*wv3));
gyR = -4.*p2.*(6.*p3.^2-1).*(-gg1.*s1.*s2.*wv2+gg2.*((s1.^2-s2.^2).*wv2-s1.*s2.*wv1)+gg3.*(s1.^2+s3.^2).*wv3)+4.*(6.*p2.^3-4.*p2).*(gg1.*(s1.^2+s2.^2).*wv1+gg3.*(s2.^2+s3.^2).*wv3+s1.*s3.*(gg3.*wv1+gg1.*wv3))+12.*p1.*p2.*p3.*(4.*gg2.*s1.*s3.*wv2-2.*s2.*s3.*(gg2.*wv1+gg1.*wv2)-2.*(s1.^2+s3.^2).*(gg3.*wv1+gg1.*wv3)-2.*s1.*s2.*(gg3.*wv2+gg2.*wv3))+2.*p1.*(3.*p2.^2-1).*(-2.*(s1.^2+s2.^2).*(gg2.*wv1+gg1.*wv2)+4.*gg3.*s1.*s2.*wv3-2.*s2.*s3.*(gg3.*wv1+gg1.*wv3)-2.*s1.*s3.*(gg3.*wv2+gg2.*wv3))+p1.*(6.*p2.^2-1).*(-2.*(s1.^2+s2.^2).*(gg2.*wv1+gg1.*wv2)+4.*gg3.*s1.*s2.*wv3-2.*s2.*s3.*(gg3.*wv1+gg1.*wv3)-2.*s1.*s3.*(gg3.*wv2+gg2.*wv3))+2.*(3.*p2.^2-1).*p3.*(4.*gg1.*s2.*s3.*wv1-2.*s1.*s3.*(gg2.*wv1+gg1.*wv2)-2.*s1.*s2.*(gg3.*wv1+gg1.*wv3)-2.*(s2.^2+s3.^2).*(gg3.*wv2+gg2.*wv3))+(6.*p2.^2-1).*p3.*(4.*gg1.*s2.*s3.*wv1-2.*s1.*s3.*(gg2.*wv1+gg1.*wv2)-2.*s1.*s2.*(gg3.*wv1+gg1.*wv3)-2.*(s2.^2+s3.^2).*(gg3.*wv2+gg2.*wv3))+4.*(6.*p1.^2-1).*p2.*(gg2.*((s1.^2+s2.^2).*wv2+s2.*s3.*wv3)+gg3.*(s2.*s3.*wv2+(s1.^2+s3.^2).*wv3));
eMFTxy = (cosd(inc)).*cos(Alpha2).*2.*alpha_0.^2./beta_0.^2.*exR+(cosd(inc)).*sin(Alpha).*2.*alpha_0.^2./beta_0.^2.*eyR;
gMFTxy = 2.*((cosd(inc)).*cos(Alpha2).*gxR+(cosd(inc)).*sin(Alpha2).*gyR);
% Far field
exR = -2.*(p1.^3.*(s1.^2-1)+2.*p1.^2.*s1.*(p2.*s2+p3.*s3)+p1.*(p2.^2.*(s2.^2-1)+2.*p3.*p2.*s2.*s3+p3.^2.*(s3.^2-1)-s1.^2+1)-...
    s1.*(p2.*s2+p3.*s3)).*(-gg1.*s1.*(s2.*wv2+s3.*wv3)+gg3.*((s1.^2-s3.^2).*wv3-s3.*(s1.*wv1+s2.*wv2))+...
    gg2.*(s1.^2.*wv2-s1.*s2.*wv1-s2.*(s2.*wv2+s3.*wv3)));
eyR = -2.*(p2.^3.*(s2.^2-1)+2.*p2.^2.*s2.*(p1.*s1+p3.*s3)+p2.*(p1.^2.*(s1.^2-1)+2.*p3.*p1.*s1.*s3+p3.^2.*(s3.^2-1)-s2.^2+1)-...
    s2.*(p1.*s1+p3.*s3)).*(-gg1.*s1.*(s2.*wv2+s3.*wv3)+gg3.*((s1.^2-s3.^2).*wv3-s3.*(s1.*wv1+s2.*wv2))+...
    gg2.*(s1.^2.*wv2-s1.*s2.*wv1-s2.*(s2.*wv2+s3.*wv3)));
gxR = p3.*p1.^2.*(4.*gg2.*s1.*s3.*wv2-2.*s2.*s3.*(gg2.*wv1+gg1.*wv2)-2.*(s1.^2+s3.^2).*(gg3.*wv1+gg1.*wv3)-...
    2.*s1.*s2.*(gg3.*wv2+gg2.*wv3))+p2.*p1.^2.*(-2.*(s1.^2+s2.^2).*(gg2.*wv1+gg1.*wv2)+4.*gg3.*s1.*s2.*wv3-...
    2.*s2.*s3.*(gg3.*wv1+gg1.*wv3)-2.*s1.*s3.*(gg3.*wv2+gg2.*wv3))+4.*p3.^2.*p1.*(gg1.*(s1.^2+s3.^2).*wv1+...
    gg2.*(s2.^2+s3.^2).*wv2+s1.*s2.*(gg2.*wv1+gg1.*wv2))+4.*p2.^2.*p1.*(gg1.*(s1.^2+s2.^2).*wv1+gg3.*(s2.^2+s3.^2).*wv3+...
    s1.*s3.*(gg3.*wv1+gg1.*wv3))+2.*p2.*p3.*p1.*(4.*gg1.*s2.*s3.*wv1-2.*s1.*s3.*(gg2.*wv1+gg1.*wv2)-2.*s1.*s2.*(gg3.*wv1+gg1.*wv3)-...
    2.*(s2.^2+s3.^2).*(gg3.*wv2+gg2.*wv3))+4.*(p1.^2-1).*p1.*(gg2.*((s1.^2+s2.^2).*wv2+s2.*s3.*wv3)+gg3.*(s2.*s3.*wv2+(s1.^2+s3.^2).*wv3))+...
    (p1.^2-1).*p3.*(4.*gg2.*s1.*s3.*wv2-2.*s2.*s3.*(gg2.*wv1+gg1.*wv2)-2.*(s1.^2+s3.^2).*(gg3.*wv1+gg1.*wv3)-2.*s1.*s2.*(gg3.*wv2+gg2.*wv3))+...
    (p1.^2-1).*p2.*(-2.*(s1.^2+s2.^2).*(gg2.*wv1+gg1.*wv2)+4.*gg3.*s1.*s2.*wv3-2.*s2.*s3.*(gg3.*wv1+gg1.*wv3)-2.*s1.*s3.*(gg3.*wv2+gg2.*wv3));
gyR = 4.*p2.*p1.^2.*(gg2.*((s1.^2+s2.^2).*wv2+s2.*s3.*wv3)+gg3.*(s2.*s3.*wv2+(s1.^2+s3.^2).*wv3))+...
    2.*p2.*p3.*p1.*(4.*gg2.*s1.*s3.*wv2-2.*s2.*s3.*(gg2.*wv1+gg1.*wv2)-2.*(s1.^2+s3.^2).*(gg3.*wv1+gg1.*wv3)-...
    2.*s1.*s2.*(gg3.*wv2+gg2.*wv3))+p2.^2.*p1.*(-2.*(s1.^2+s2.^2).*(gg2.*wv1+gg1.*wv2)+4.*gg3.*s1.*s2.*wv3-...
    2.*s2.*s3.*(gg3.*wv1+gg1.*wv3)-2.*s1.*s3.*(gg3.*wv2+gg2.*wv3))+(p2.^2-1).*p1.*(-2.*(s1.^2+s2.^2).*(gg2.*wv1+gg1.*wv2)+...
    4.*gg3.*s1.*s2.*wv3-2.*s2.*s3.*(gg3.*wv1+gg1.*wv3)-2.*s1.*s3.*(gg3.*wv2+gg2.*wv3))+4.*p2.*p3.^2.*(gg1.*(s1.^2+s3.^2).*wv1+...
    gg2.*(s2.^2+s3.^2).*wv2+s1.*s2.*(gg2.*wv1+gg1.*wv2))+4.*p2.*(p2.^2-1).*(gg1.*(s1.^2+s2.^2).*wv1+gg3.*(s2.^2+s3.^2).*wv3+...
    s1.*s3.*(gg3.*wv1+gg1.*wv3))+p2.^2.*p3.*(4.*gg1.*s2.*s3.*wv1-2.*s1.*s3.*(gg2.*wv1+gg1.*wv2)-2.*s1.*s2.*(gg3.*wv1+gg1.*wv3)-...
    2.*(s2.^2+s3.^2).*(gg3.*wv2+gg2.*wv3))+(p2.^2-1).*p3.*(4.*gg1.*s2.*s3.*wv1-2.*s1.*s3.*(gg2.*wv1+gg1.*wv2)-2.*s1.*s2.*(gg3.*wv1+gg1.*wv3)-...
    2.*(s2.^2+s3.^2).*(gg3.*wv2+gg2.*wv3));
eFFTxy = (cosd(inc)).*cos(Alpha2).*2.*alpha_0.^2./beta_0.^2.*exR+(cosd(inc)).*sin(Alpha).*2.*alpha_0.^2./beta_0.^2.*eyR;
gFFTxy = 2.*((cosd(inc)).*cos(Alpha2).*gxR+(cosd(inc)).*sin(Alpha2).*gyR);

%The elasticity tensor for a single crystal olivine 
rho=3300;
a = 0.21;
d_a = a.*alpha_0;
d_b = 1.25.*a.*beta_0;
aa = rho.*(alpha_0 - d_a./2).^2;
cc = rho.*(alpha_0 + d_a./2).^2;
ll = rho.*(beta_0 + d_b./2).^2;
nn = rho.*(beta_0 - d_b./2).^2;
ac = rho.*(alpha_0).^2;
ff = -ll + sqrt((2.*ac).^2 - 2.*ac.*(aa + cc + 2.*ll) + (aa + ll).*(cc + ll));

C33 = aa./(10.^9);
C11 = cc./(10.^9);
C13 = ff./(10.^9);
C44 = nn./(10.^9);
C66 = ll./(10.^9);
% E. H. Abramson, J. M. Brown, L. J. Slutsky, J. Zaug; 1997
% C11=320.5;
% C12=68.1;
% C13=71.6;
% C22=196.5;
% C23=76.8;
% C33=233.5;
% C44=64.0;
% C55=77.0;
% C66=78.7;
% Mondal & Long; 2019
% C11=173;
% C12=69;
% C13=58;
% C33=272;
% C44=61;
% C66=52;
eps = (10.^9.*(C11-C33)./(2.*rho.*alpha_0.^2));
del = (10.^9.*(C13-C33+2*C44)./(rho.*alpha_0.^2));
gam = (10.^9.*(C66-C44)./(2.*rho.*beta_0.^2)); 
dFFTxy = -eFFTxy;
dMFTxy = -eMFTxy;
dNFTxy = -eNFTxy;
gam = -gam;
eps = -eps;

ggTr = (eps.*(eNFTxy.*gNF+eMFTxy.*gMF+eFFTxy.*gFF)+del.*(dNFTxy.*gNF+dMFTxy.*gMF+dFFTxy.*gFF)+gam.*(gNFTxy.*gNF+gMFTxy.*gMF+gFFTxy.*gFF)).*strength_in;

% calculate splitting intensity (correlation)
s=ggTr.*dx.*dy.*dz;
si = sum(s(:));
vi = s;

end
