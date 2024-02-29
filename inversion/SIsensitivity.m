function [WF,EKx,EKphi,EKtheta,EKx2,EKphi2,EKtheta2,timestep]=SIsensitivity(data,indoi,w,sflag,vflag,dec,strength_in,phi_in,theta_in,zlim,nx,ny,nz,edgeadd,nny)
% Copyright 2024 F.Link and M.D.Long 

if nargin<13
    edgeadd = 100000;
    nny = 4;
end

AZIMUTH=data.baz(indoi);
per=data.per(indoi);
par=data.p(indoi);

X_ST=data.x.*1000;
Y_ST=data.y.*1000;

load('ak135.mat','speed');
depth=speed(:,1).*1000; %convert everything to SI units
vp=speed(:,3).*1000;
vs=speed(:,2).*1000;
NN=max(size(vp(:))); 
datarep=[ones(NN,1) depth(:)];
coeff_p=pinv(datarep'*datarep)*datarep'*vp(:);
coeff_s=pinv(datarep'*datarep)*datarep'*vs(:);

% allocate sensitivity kernels
EKx = zeros(nx,ny,nz);
EKphi = EKx;
EKtheta = EKx;
EKx2 = EKx;
EKphi2 = EKphi;
EKtheta2 = EKtheta;

X=abs(max(X_ST)-min(X_ST)+edgeadd);
dx=X/nx;
if ny == 1
    Y_ST = zeros(size(X_ST));
    dy = nny*dx;
    y = -(edgeadd/2):dy:(edgeadd/2);
    strength_in = repmat(strength_in,1,length(y),1);
    phi_in = repmat(phi_in,1,length(y),1);
    theta_in = repmat(theta_in,1,length(y),1);
else  
    y=linspace(min(Y_ST)-(edgeadd/2),max(Y_ST)+(edgeadd/2),ny);
end
x=linspace(min(X_ST)-(edgeadd/2),max(X_ST)+(edgeadd/2),nx);
z=linspace(zlim(1),zlim(2),nz);
X_ST = X_ST(indoi);
Y_ST = Y_ST(indoi);

M1=max(size(x(:)));
M2=max(size(y(:)));
M3=max(size(z(:)));

vpe=coeff_p(1)+coeff_p(2).*z; %smoothed out linear P-wave speed
vse=coeff_s(1)+coeff_s(2).*z; %smoothed out linear S-wave speed and so on

%speed profile extracted from the ak135 model
alpha_0 = ones(M1,M2,1).*reshape(vpe,1,1,[]);
beta_0 = ones(M1,M2,1).*reshape(vse,1,1,[]);
if M2 == 1 || M1 == 1
    alpha_0 = reshape(alpha_0,M1,M2,M3);
    beta_0 = reshape(beta_0,M1,M2,M3);
end

% initialize wave
F_max = 1;
F_min=0.01;
df=(F_max-F_min)/80;
f=F_min:df:F_max;
omega=2.*pi.*f;

% loop over events
Nst = length(X_ST);

% position vector
[Y,X,Z] = meshgrid(y,x,z);

% initialize variable
for i = 1:Nst
    WF(i).vals = 0;
    WF(i).misfit = 0;
    DT(i).vals = data.si(indoi(i));
end

timeref = tic;
if isempty(gcp)
    parpool(feature('NumCores'));
end
% parfor iiii = 1:Nst
parfor iiii = 1:Nst

% station location
x_1=X_ST(iiii);
y_1=Y_ST(iiii);
z_1=0;
% keyboard
[si,dSKxTr,dSKphiTr,dSKthetaTr] = ForwardSIsensitivity(per(iiii),par(iiii),AZIMUTH(iiii),x_1,y_1,z_1,X,Y,Z,omega,alpha_0,beta_0,strength_in,phi_in,theta_in,dec,vflag);

WF(iiii).vals = si;

if sflag
    WF(iiii).misfit = 1;
else
    D_in = DT(iiii).vals;
    WF(iiii).misfit = si-D_in;
end

SKx = dSKxTr.*WF(iiii).misfit;
SKphi = dSKphiTr.*WF(iiii).misfit;
SKtheta = dSKthetaTr.*WF(iiii).misfit;

if ny == 1
    temp1 = reshape(mean(real(SKx),2),M1,1,M3);
    temp2 = reshape(mean(real(SKphi),2),M1,1,M3);
    temp3 = reshape(mean(real(SKtheta),2),M1,1,M3);
else
    temp1 = real(SKx);
    temp2 = real(SKphi);
    temp3 = real(SKtheta);
end
EKx = EKx+temp1.*w(iiii);
EKphi = EKphi+temp2.*w(iiii);
EKtheta = EKtheta+temp3.*w(iiii);

temp = zeros(size(temp1));
temp(max(abs(temp1)).*0.1<abs(temp1)) = 1;
EKx2 = EKx2+temp.*w(iiii);
temp = zeros(size(temp2));
temp(max(abs(temp2)).*0.1<abs(temp2)) = 1;
EKphi2 = EKphi2+temp.*w(iiii);
temp = zeros(size(temp3));
temp(max(abs(temp3)).*0.1<abs(temp3)) = 1;
EKtheta2 = EKtheta2+temp.*w(iiii);
end
timestep = toc(timeref);


end

