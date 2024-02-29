function [dSKthetaTr] = ForwardSIsensitivityTheta(per,par,AZIMUTH,x_1,y_1,z_1,X,Y,Z,omega,alpha_0,beta_0,strength_in,phi_in,theta_in,vflag)
% Copyright 2024 F.Link and M.D.Long 

MM = length(X);

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
wv(2,:) = sind(inc(:)).*cosd(AZIMUTH+180);
wv(3,:) = cosd(inc(:));

% detour of scattered wave relative to direct wave
D = abs(dot(wv,[X-x_1 Y-y_1 Z-z_1]',1)');

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
M = iota.*exp((iota.*(R-D)./beta_0).*reshape(omega,1,[])).*(ones(MM,1).*reshape(m,1,[])).*dom;
sumsum = sum(omega.^2.*m.*dom);
fFF = reshape(omega.^3,1,[])./4./pi./beta_0.^2./R;
fMFF = sum(real(fFF.*M),2);
gFF=real(fMFF)./(sumsum);

fMF = -iota.*reshape(omega.^2,1,[])./4./pi./beta_0./R.^2;
fMMF = sum(real(fMF.*M),2);
gMF=real(fMMF)./(sumsum);

fNF = -1./4./pi./R.^3;
fMNF = fNF.*sum(real(M).*reshape(omega,1,[]),2);
gNF=real(fMNF)./(sumsum);

% Radiation pattern
wv1 = wv(1,:)';
wv2 = wv(2,:)';
wv3 = wv(3,:)';
gg1 = gg(1,:)';
gg2 = gg(2,:)';
gg3 = gg(3,:)';

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
eps = (10.^9.*(C11-C33)./(2.*rho.*alpha_0.^2));
del = (10.^9.*(C13-C33+2*C44)./(rho.*alpha_0.^2));
gam = (10.^9.*(C66-C44)./(2.*rho.*beta_0.^2)); 
eps = -eps;
gam = -gam;

% Prepare sensitivity kernels for plunge of anisotropic tensor
% derivative for theta > s1 = s1b; s2 = s2b; s3 = s3b >> d(...)/d_theta
s1b = sin(-theta_in).*sin(phi_in);
s2b = -sin(-theta_in).*cos(phi_in);
s3b = cos(-theta_in);

% Near Field
exR = 6.*(-5.*p1.*(p1.^2+p2.^2+p3.^2-1)+2.*(5.*p1.^2-1).*s1.*(p2.*s2+p3.*s3)+p1.*(5.*p1.^2-3).*s1.^2+p1.*((5.*p2.^2-1).*s2.^2+10.*p2.*p3.*s2.*s3+(5.*p3.^2-1).*s3.^2)).*(s1.*(s2b.*(gg1.*wv2+gg2.*wv1)+s3b.*(gg1.*wv3+gg3.*wv1))+s1b.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))-2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+(gg2.*s2b+gg3.*s3b).*(wv2.*s2+wv3.*s3)+(gg2.*s2+gg3.*s3).*(wv2.*s2b+wv3.*s3b))+6.*(s1.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))-(s1.^2.*(gg2.*wv2+gg3.*wv3))+(gg2.*s2+gg3.*s3).*(wv2.*s2+wv3.*s3)).*(2.*(5.*p1.^2-1).*s1.*(p2.*s2b+p3.*s3b)+2.*(5.*p1.^2-1).*s1b.*(p2.*s2+p3.*s3)+2.*p1.*(5.*p1.^2-3).*s1.*s1b+p1.*(2.*(5.*p2.^2-1).*s2.*s2b+10.*p2.*p3.*s2.*s3b+10.*p2.*p3.*s2b.*s3+2.*(5.*p3.^2-1).*s3.*s3b));
eyR = 6.*(s1.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))-(s1.^2.*(gg2.*wv2+gg3.*wv3))+(gg2.*s2+gg3.*s3).*(wv2.*s2+wv3.*s3)).*(2.*(5.*p1.^2-1).*p2.*s1.*s1b+2.*p1.*s1.*((5.*p2.^2-1).*s2b+5.*p2.*p3.*s3b)+2.*p1.*s1b.*((5.*p2.^2-1).*s2+5.*p2.*p3.*s3)+2.*(5.*p2.^2-1).*p3.*s2.*s3b+2.*(5.*p2.^2-1).*p3.*s2b.*s3+2.*p2.*(5.*p2.^2-3).*s2.*s2b+2.*p2.*(5.*p3.^2-1).*s3.*s3b)+6.*(-5.*p2.*(p1.^2+p2.^2+p3.^2-1)+(5.*p1.^2-1).*p2.*s1.^2+2.*p1.*s1.*((5.*p2.^2-1).*s2+5.*p2.*p3.*s3)+2.*(5.*p2.^2-1).*p3.*s2.*s3+p2.*(5.*p2.^2-3).*s2.^2+p2.*(5.*p3.^2-1).*s3.^2).*(s1.*(s2b.*(gg1.*wv2+gg2.*wv1)+s3b.*(gg1.*wv3+gg3.*wv1))+s1b.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))-2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+(gg2.*s2b+gg3.*s3b).*(wv2.*s2+wv3.*s3)+(gg2.*s2+gg3.*s3).*(wv2.*s2b+wv3.*s3b));
gxR = 6.*(-2.*(5.*p1.^2-1).*p2.*(s2b.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))+s2.*(s2b.*(gg1.*wv2+gg2.*wv1)+s3b.*(gg1.*wv3+gg3.*wv1))+2.*s1.*s1b.*(gg1.*wv2+gg2.*wv1)+s1.*(s3b.*(gg2.*wv3+gg3.*wv2)-2.*gg3.*wv3.*s2b)+s1b.*(s3.*(gg2.*wv3+gg3.*wv2)-2.*gg3.*wv3.*s2))-2.*(5.*p1.^2-1).*p3.*(s3b.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))+s3.*(s2b.*(gg1.*wv2+gg2.*wv1)+s3b.*(gg1.*wv3+gg3.*wv1))+2.*s1.*s1b.*(gg1.*wv3+gg3.*wv1)+s1.*(s2b.*(gg2.*wv3+gg3.*wv2)-2.*gg2.*wv2.*s3b)+s1b.*(s2.*(gg2.*wv3+gg3.*wv2)-2.*gg2.*wv2.*s3))+2.*p1.*(5.*p2.^2-1).*((gg1.*s1b+gg3.*s3b).*(wv1.*s1+wv3.*s3)+(gg1.*s1+gg3.*s3).*(wv1.*s1b+wv3.*s3b)-2.*gg2.*wv2.*s2.*s2b)+5.*p1.*p2.*p3.*(-2.*s1.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s1b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*s1.*s2b.*(gg1.*wv3+gg3.*wv1)-2.*s1b.*s2.*(gg1.*wv3+gg3.*wv1)+4.*gg1.*wv1.*s2.*s3b+4.*gg1.*wv1.*s2b.*s3-2.*(gg2.*wv3+gg3.*wv2).*(2.*s2.*s2b+2.*s3.*s3b))+2.*p1.*(5.*p3.^2-1).*(s1.*s2b.*(gg1.*wv2+gg2.*wv1)+s1b.*s2.*(gg1.*wv2+gg2.*wv1)-2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+2.*gg2.*wv2.*s2.*s2b-2.*gg3.*wv3.*s3.*s3b)+2.*p1.*(5.*p1.^2-3).*(2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+(gg2.*s2b+gg3.*s3b).*(wv2.*s2+wv3.*s3)+(gg2.*s2+gg3.*s3).*(wv2.*s2b+wv3.*s3b)));
gyR = 6.*(-2.*p1.*(5.*p2.^2-1).*(s2b.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))+s2.*(s2b.*(gg1.*wv2+gg2.*wv1)+s3b.*(gg1.*wv3+gg3.*wv1))+2.*s1.*s1b.*(gg1.*wv2+gg2.*wv1)+s1.*(s3b.*(gg2.*wv3+gg3.*wv2)-2.*gg3.*wv3.*s2b)+s1b.*(s3.*(gg2.*wv3+gg3.*wv2)-2.*gg3.*wv3.*s2))+5.*p1.*p2.*p3.*(-2.*s2.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s2b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*(gg1.*wv3+gg3.*wv1).*(2.*s1.*s1b+2.*s3.*s3b)-2.*s1.*s2b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s2.*(gg2.*wv3+gg3.*wv2)+4.*gg2.*wv2.*s1.*s3b+4.*gg2.*wv2.*s1b.*s3)+(5.*p2.^2-1).*p3.*(-2.*s1.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s1b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*s1.*s2b.*(gg1.*wv3+gg3.*wv1)-2.*s1b.*s2.*(gg1.*wv3+gg3.*wv1)+4.*gg1.*wv1.*s2.*s3b+4.*gg1.*wv1.*s2b.*s3-2.*(gg2.*wv3+gg3.*wv2).*(2.*s2.*s2b+2.*s3.*s3b))+2.*p2.*(5.*p2.^2-3).*((gg1.*s1b+gg3.*s3b).*(wv1.*s1+wv3.*s3)+(gg1.*s1+gg3.*s3).*(wv1.*s1b+wv3.*s3b)-2.*gg2.*wv2.*s2.*s2b)+2.*p2.*(5.*p3.^2-1).*(s1.*s2b.*(gg1.*wv2+gg2.*wv1)+s1b.*s2.*(gg1.*wv2+gg2.*wv1)-2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+2.*gg2.*wv2.*s2.*s2b-2.*gg3.*wv3.*s3.*s3b)+2.*(5.*p1.^2-1).*p2.*(2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+(gg2.*s2b+gg3.*s3b).*(wv2.*s2+wv3.*s3)+(gg2.*s2+gg3.*s3).*(wv2.*s2b+wv3.*s3b)));
eNFTxy = (cosd(inc)).*cos(Alpha2).*2.*alpha_0.^2./beta_0.^2.*exR+(cosd(inc)).*sin(Alpha).*2.*alpha_0.^2./beta_0.^2.*eyR;
gNFTxy = 2.*((cosd(inc)).*cos(Alpha2).*gxR+(cosd(inc)).*sin(Alpha2).*gyR);
% Mid field
exR = 2.*((6.*p1.^3-4.*p1).*s1.^2-6.*p1.*(p1.^2+p2.^2+p3.^2-1)+3.*(4.*p1.^2-1).*s1.*(p2.*s2+p3.*s3)+p1.*((6.*p2.^2-1).*s2.^2+12.*p2.*p3.*s2.*s3+(6.*p3.^2-1).*s3.^2)).*(s1.*(s2b.*(gg1.*wv2+gg2.*wv1)+s3b.*(gg1.*wv3+gg3.*wv1))+s1b.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))-2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+(gg2.*s2b+gg3.*s3b).*(wv2.*s2+wv3.*s3)+(gg2.*s2+gg3.*s3).*(wv2.*s2b+wv3.*s3b))+2.*(s1.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))-(s1.^2.*(gg2.*wv2+gg3.*wv3))+(gg2.*s2+gg3.*s3).*(wv2.*s2+wv3.*s3)).*(2.*(6.*p1.^3-4.*p1).*s1.*s1b+3.*(4.*p1.^2-1).*s1.*(p2.*s2b+p3.*s3b)+3.*(4.*p1.^2-1).*s1b.*(p2.*s2+p3.*s3)+p1.*(2.*(6.*p2.^2-1).*s2.*s2b+12.*p2.*p3.*s2.*s3b+12.*p2.*p3.*s2b.*s3+2.*(6.*p3.^2-1).*s3.*s3b));
eyR = 2.*(s1.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))-(s1.^2.*(gg2.*wv2+gg3.*wv3))+(gg2.*s2+gg3.*s3).*(wv2.*s2+wv3.*s3)).*(2.*(6.*p1.^2-1).*p2.*s1.*s1b+3.*p1.*s1.*((4.*p2.^2-1).*s2b+4.*p2.*p3.*s3b)+3.*p1.*s1b.*((4.*p2.^2-1).*s2+4.*p2.*p3.*s3)+2.*(6.*p2.^3-4.*p2).*s2.*s2b+3.*(4.*p2.^2-1).*p3.*s2.*s3b+3.*(4.*p2.^2-1).*p3.*s2b.*s3+2.*p2.*(6.*p3.^2-1).*s3.*s3b)+2.*(-6.*p2.*(p1.^2+p2.^2+p3.^2-1)+(6.*p1.^2-1).*p2.*s1.^2+3.*p1.*s1.*((4.*p2.^2-1).*s2+4.*p2.*p3.*s3)+(6.*p2.^3-4.*p2).*s2.^2+3.*(4.*p2.^2-1).*p3.*s2.*s3+p2.*(6.*p3.^2-1).*s3.^2).*(s1.*(s2b.*(gg1.*wv2+gg2.*wv1)+s3b.*(gg1.*wv3+gg3.*wv1))+s1b.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))-2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+(gg2.*s2b+gg3.*s3b).*(wv2.*s2+wv3.*s3)+(gg2.*s2+gg3.*s3).*(wv2.*s2b+wv3.*s3b));
gxR = 2.*(3.*p1.^2-1).*p2.*(-2.*(gg1.*wv2+gg2.*wv1).*(2.*s1.*s1b+2.*s2.*s2b)-2.*s2.*s3b.*(gg1.*wv3+gg3.*wv1)-2.*s2b.*s3.*(gg1.*wv3+gg3.*wv1)-2.*s1.*s3b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s3.*(gg2.*wv3+gg3.*wv2)+4.*gg3.*wv3.*s1.*s2b+4.*gg3.*wv3.*s1b.*s2)+(6.*p1.^2-1).*p2.*(-2.*(gg1.*wv2+gg2.*wv1).*(2.*s1.*s1b+2.*s2.*s2b)-2.*s2.*s3b.*(gg1.*wv3+gg3.*wv1)-2.*s2b.*s3.*(gg1.*wv3+gg3.*wv1)-2.*s1.*s3b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s3.*(gg2.*wv3+gg3.*wv2)+4.*gg3.*wv3.*s1.*s2b+4.*gg3.*wv3.*s1b.*s2)+2.*(3.*p1.^2-1).*p3.*(-2.*s2.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s2b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*(gg1.*wv3+gg3.*wv1).*(2.*s1.*s1b+2.*s3.*s3b)-2.*s1.*s2b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s2.*(gg2.*wv3+gg3.*wv2)+4.*gg2.*wv2.*s1.*s3b+4.*gg2.*wv2.*s1b.*s3)+(6.*p1.^2-1).*p3.*(-2.*s2.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s2b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*(gg1.*wv3+gg3.*wv1).*(2.*s1.*s1b+2.*s3.*s3b)-2.*s1.*s2b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s2.*(gg2.*wv3+gg3.*wv2)+4.*gg2.*wv2.*s1.*s3b+4.*gg2.*wv2.*s1b.*s3)+4.*p1.*(6.*p2.^2-1).*((gg1.*s1b+gg3.*s3b).*(wv1.*s1+wv3.*s3)+(gg1.*s1+gg3.*s3).*(wv1.*s1b+wv3.*s3b)-2.*gg2.*wv2.*s2.*s2b)+12.*p1.*p2.*p3.*(-2.*s1.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s1b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*s1.*s2b.*(gg1.*wv3+gg3.*wv1)-2.*s1b.*s2.*(gg1.*wv3+gg3.*wv1)+4.*gg1.*wv1.*s2.*s3b+4.*gg1.*wv1.*s2b.*s3-2.*(gg2.*wv3+gg3.*wv2).*(2.*s2.*s2b+2.*s3.*s3b))+4.*p1.*(6.*p3.^2-1).*(s1.*s2b.*(gg1.*wv2+gg2.*wv1)+s1b.*s2.*(gg1.*wv2+gg2.*wv1)-2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+2.*gg2.*wv2.*s2.*s2b-2.*gg3.*wv3.*s3.*s3b)+8.*p1.*(3.*p1.^2-2).*(2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+(gg2.*s2b+gg3.*s3b).*(wv2.*s2+wv3.*s3)+(gg2.*s2+gg3.*s3).*(wv2.*s2b+wv3.*s3b));
gyR = 2.*p1.*(3.*p2.^2-1).*(-2.*(gg1.*wv2+gg2.*wv1).*(2.*s1.*s1b+2.*s2.*s2b)-2.*s2.*s3b.*(gg1.*wv3+gg3.*wv1)-2.*s2b.*s3.*(gg1.*wv3+gg3.*wv1)-2.*s1.*s3b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s3.*(gg2.*wv3+gg3.*wv2)+4.*gg3.*wv3.*s1.*s2b+4.*gg3.*wv3.*s1b.*s2)+p1.*(6.*p2.^2-1).*(-2.*(gg1.*wv2+gg2.*wv1).*(2.*s1.*s1b+2.*s2.*s2b)-2.*s2.*s3b.*(gg1.*wv3+gg3.*wv1)-2.*s2b.*s3.*(gg1.*wv3+gg3.*wv1)-2.*s1.*s3b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s3.*(gg2.*wv3+gg3.*wv2)+4.*gg3.*wv3.*s1.*s2b+4.*gg3.*wv3.*s1b.*s2)+12.*p1.*p2.*p3.*(-2.*s2.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s2b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*(gg1.*wv3+gg3.*wv1).*(2.*s1.*s1b+2.*s3.*s3b)-2.*s1.*s2b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s2.*(gg2.*wv3+gg3.*wv2)+4.*gg2.*wv2.*s1.*s3b+4.*gg2.*wv2.*s1b.*s3)+2.*(3.*p2.^2-1).*p3.*(-2.*s1.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s1b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*s1.*s2b.*(gg1.*wv3+gg3.*wv1)-2.*s1b.*s2.*(gg1.*wv3+gg3.*wv1)+4.*gg1.*wv1.*s2.*s3b+4.*gg1.*wv1.*s2b.*s3-2.*(gg2.*wv3+gg3.*wv2).*(2.*s2.*s2b+2.*s3.*s3b))+(6.*p2.^2-1).*p3.*(-2.*s1.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s1b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*s1.*s2b.*(gg1.*wv3+gg3.*wv1)-2.*s1b.*s2.*(gg1.*wv3+gg3.*wv1)+4.*gg1.*wv1.*s2.*s3b+4.*gg1.*wv1.*s2b.*s3-2.*(gg2.*wv3+gg3.*wv2).*(2.*s2.*s2b+2.*s3.*s3b))+8.*p2.*(3.*p2.^2-2).*((gg1.*s1b+gg3.*s3b).*(wv1.*s1+wv3.*s3)+(gg1.*s1+gg3.*s3).*(wv1.*s1b+wv3.*s3b)-2.*gg2.*wv2.*s2.*s2b)+4.*p2.*(6.*p3.^2-1).*(s1.*s2b.*(gg1.*wv2+gg2.*wv1)+s1b.*s2.*(gg1.*wv2+gg2.*wv1)-2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+2.*gg2.*wv2.*s2.*s2b-2.*gg3.*wv3.*s3.*s3b)+4.*(6.*p1.^2-1).*p2.*(2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+(gg2.*s2b+gg3.*s3b).*(wv2.*s2+wv3.*s3)+(gg2.*s2+gg3.*s3).*(wv2.*s2b+wv3.*s3b));
eMFTxy = (cosd(inc)).*cos(Alpha2).*2.*alpha_0.^2./beta_0.^2.*exR+(cosd(inc)).*sin(Alpha).*2.*alpha_0.^2./beta_0.^2.*eyR;
gMFTxy = 2.*((cosd(inc)).*cos(Alpha2).*gxR+(cosd(inc)).*sin(Alpha2).*gyR);
% Far field
exR = 2.*((p1.*s1+p2.*s2+p3.*s3).*((p1.^2-1).*s1+p1.*(p2.*s2+p3.*s3))-p1.*(p1.^2+p2.^2+p3.^2-1)).*(s1.*(s2b.*(gg1.*wv2+gg2.*wv1)+s3b.*(gg1.*wv3+gg3.*wv1))+s1b.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))-2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+(gg2.*s2b+gg3.*s3b).*(wv2.*s2+wv3.*s3)+(gg2.*s2+gg3.*s3).*(wv2.*s2b+wv3.*s3b))+2.*(s1.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))-(s1.^2.*(gg2.*wv2+gg3.*wv3))+(gg2.*s2+gg3.*s3).*(wv2.*s2+wv3.*s3)).*(((p1.^2-1).*s1+p1.*(p2.*s2+p3.*s3)).*(p1.*s1b+p2.*s2b+p3.*s3b)+(p1.*s1+p2.*s2+p3.*s3).*((p1.^2-1).*s1b+p1.*(p2.*s2b+p3.*s3b)));
eyR = 2.*((p1.*s1+p2.*s2+p3.*s3).*(p1.*p2.*s1+(p2.^2-1).*s2+p2.*p3.*s3)-p2.*(p1.^2+p2.^2+p3.^2-1)).*(s1.*(s2b.*(gg1.*wv2+gg2.*wv1)+s3b.*(gg1.*wv3+gg3.*wv1))+s1b.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))-2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+(gg2.*s2b+gg3.*s3b).*(wv2.*s2+wv3.*s3)+(gg2.*s2+gg3.*s3).*(wv2.*s2b+wv3.*s3b))+2.*(s1.*(s2.*(gg1.*wv2+gg2.*wv1)+s3.*(gg1.*wv3+gg3.*wv1))-(s1.^2.*(gg2.*wv2+gg3.*wv3))+(gg2.*s2+gg3.*s3).*(wv2.*s2+wv3.*s3)).*((p1.*p2.*s1+(p2.^2-1).*s2+p2.*p3.*s3).*(p1.*s1b+p2.*s2b+p3.*s3b)+(p1.*s1+p2.*s2+p3.*s3).*(p1.*p2.*s1b+(p2.^2-1).*s2b+p2.*p3.*s3b));
gxR = p1.^2.*p2.*(-2.*(gg1.*wv2+gg2.*wv1).*(2.*s1.*s1b+2.*s2.*s2b)-2.*s2.*s3b.*(gg1.*wv3+gg3.*wv1)-2.*s2b.*s3.*(gg1.*wv3+gg3.*wv1)-2.*s1.*s3b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s3.*(gg2.*wv3+gg3.*wv2)+4.*gg3.*wv3.*s1.*s2b+4.*gg3.*wv3.*s1b.*s2)+(p1.^2-1).*p2.*(-2.*(gg1.*wv2+gg2.*wv1).*(2.*s1.*s1b+2.*s2.*s2b)-2.*s2.*s3b.*(gg1.*wv3+gg3.*wv1)-2.*s2b.*s3.*(gg1.*wv3+gg3.*wv1)-2.*s1.*s3b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s3.*(gg2.*wv3+gg3.*wv2)+4.*gg3.*wv3.*s1.*s2b+4.*gg3.*wv3.*s1b.*s2)+p1.^2.*p3.*(-2.*s2.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s2b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*(gg1.*wv3+gg3.*wv1).*(2.*s1.*s1b+2.*s3.*s3b)-2.*s1.*s2b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s2.*(gg2.*wv3+gg3.*wv2)+4.*gg2.*wv2.*s1.*s3b+4.*gg2.*wv2.*s1b.*s3)+(p1.^2-1).*p3.*(-2.*s2.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s2b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*(gg1.*wv3+gg3.*wv1).*(2.*s1.*s1b+2.*s3.*s3b)-2.*s1.*s2b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s2.*(gg2.*wv3+gg3.*wv2)+4.*gg2.*wv2.*s1.*s3b+4.*gg2.*wv2.*s1b.*s3)+4.*p1.*p2.^2.*((gg1.*s1b+gg3.*s3b).*(wv1.*s1+wv3.*s3)+(gg1.*s1+gg3.*s3).*(wv1.*s1b+wv3.*s3b)-2.*gg2.*wv2.*s2.*s2b)+2.*p1.*p2.*p3.*(-2.*s1.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s1b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*s1.*s2b.*(gg1.*wv3+gg3.*wv1)-2.*s1b.*s2.*(gg1.*wv3+gg3.*wv1)+4.*gg1.*wv1.*s2.*s3b+4.*gg1.*wv1.*s2b.*s3-2.*(gg2.*wv3+gg3.*wv2).*(2.*s2.*s2b+2.*s3.*s3b))+4.*p1.*p3.^2.*((gg1.*s1b+gg2.*s2b).*(wv1.*s1+wv2.*s2)+(gg1.*s1+gg2.*s2).*(wv1.*s1b+wv2.*s2b)-2.*gg3.*wv3.*s3.*s3b)+4.*(p1.^2-1).*p1.*(2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+(gg2.*s2b+gg3.*s3b).*(wv2.*s2+wv3.*s3)+(gg2.*s2+gg3.*s3).*(wv2.*s2b+wv3.*s3b));
gyR = p1.*p2.^2.*(-2.*(gg1.*wv2+gg2.*wv1).*(2.*s1.*s1b+2.*s2.*s2b)-2.*s2.*s3b.*(gg1.*wv3+gg3.*wv1)-2.*s2b.*s3.*(gg1.*wv3+gg3.*wv1)-2.*s1.*s3b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s3.*(gg2.*wv3+gg3.*wv2)+4.*gg3.*wv3.*s1.*s2b+4.*gg3.*wv3.*s1b.*s2)+p1.*(p2.^2-1).*(-2.*(gg1.*wv2+gg2.*wv1).*(2.*s1.*s1b+2.*s2.*s2b)-2.*s2.*s3b.*(gg1.*wv3+gg3.*wv1)-2.*s2b.*s3.*(gg1.*wv3+gg3.*wv1)-2.*s1.*s3b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s3.*(gg2.*wv3+gg3.*wv2)+4.*gg3.*wv3.*s1.*s2b+4.*gg3.*wv3.*s1b.*s2)+2.*p1.*p2.*p3.*(-2.*s2.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s2b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*(gg1.*wv3+gg3.*wv1).*(2.*s1.*s1b+2.*s3.*s3b)-2.*s1.*s2b.*(gg2.*wv3+gg3.*wv2)-2.*s1b.*s2.*(gg2.*wv3+gg3.*wv2)+4.*gg2.*wv2.*s1.*s3b+4.*gg2.*wv2.*s1b.*s3)+p2.^2.*p3.*(-2.*s1.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s1b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*s1.*s2b.*(gg1.*wv3+gg3.*wv1)-2.*s1b.*s2.*(gg1.*wv3+gg3.*wv1)+4.*gg1.*wv1.*s2.*s3b+4.*gg1.*wv1.*s2b.*s3-2.*(gg2.*wv3+gg3.*wv2).*(2.*s2.*s2b+2.*s3.*s3b))+(p2.^2-1).*p3.*(-2.*s1.*s3b.*(gg1.*wv2+gg2.*wv1)-2.*s1b.*s3.*(gg1.*wv2+gg2.*wv1)-2.*s1.*s2b.*(gg1.*wv3+gg3.*wv1)-2.*s1b.*s2.*(gg1.*wv3+gg3.*wv1)+4.*gg1.*wv1.*s2.*s3b+4.*gg1.*wv1.*s2b.*s3-2.*(gg2.*wv3+gg3.*wv2).*(2.*s2.*s2b+2.*s3.*s3b))+4.*p2.*(p2.^2-1).*((gg1.*s1b+gg3.*s3b).*(wv1.*s1+wv3.*s3)+(gg1.*s1+gg3.*s3).*(wv1.*s1b+wv3.*s3b)-2.*gg2.*wv2.*s2.*s2b)+4.*p2.*p3.^2.*((gg1.*s1b+gg2.*s2b).*(wv1.*s1+wv2.*s2)+(gg1.*s1+gg2.*s2).*(wv1.*s1b+wv2.*s2b)-2.*gg3.*wv3.*s3.*s3b)+4.*p1.^2.*p2.*(2.*s1.*s1b.*(gg2.*wv2+gg3.*wv3)+(gg2.*s2b+gg3.*s3b).*(wv2.*s2+wv3.*s3)+(gg2.*s2+gg3.*s3).*(wv2.*s2b+wv3.*s3b));
eFFTxy = (cosd(inc)).*cos(Alpha2).*2.*alpha_0.^2./beta_0.^2.*exR+(cosd(inc)).*sin(Alpha).*2.*alpha_0.^2./beta_0.^2.*eyR;
gFFTxy = 2.*((cosd(inc)).*cos(Alpha2).*gxR+(cosd(inc)).*sin(Alpha2).*gyR);

dFFTxy = -eFFTxy;
dMFTxy = -eMFTxy;
dNFTxy = -eNFTxy;

dSKthetaTr = (eps.*(eNFTxy.*gNF+eMFTxy.*gMF+eFFTxy.*gFF)+del.*(dNFTxy.*gNF+dMFTxy.*gMF+dFFTxy.*gFF)+gam.*(gNFTxy.*gNF+gMFTxy.*gMF+gFFTxy.*gFF));%.*R.^0.5;

end
