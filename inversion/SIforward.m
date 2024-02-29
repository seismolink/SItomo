function [WF,timestep]=SIforward(data,indoi,sflag,vflag,strength_in,phi_in,theta_in,zlim,nx,ny,nz,edgeadd,nny)
% Copyright 2024 F.Link and M.D.Long 

if nargin<12
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
    if ~sflag
        DT(i).vals = data.si(indoi(i));
    else
        DT(i).vals = 0;
    end
end
timeref = tic;
if isempty(gcp)
    parpool(feature('NumCores'));
end
parfor iiii = 1:Nst

% station location
x_1=X_ST(iiii);
y_1=Y_ST(iiii);
z_1=0;

if sum(abs(strength_in(:)))==0
    si = 0;
    vi = zeros(M1,M2,M3);
else
    [si,vi] = ForwardSI(per(iiii),par(iiii),AZIMUTH(iiii),x_1,y_1,z_1,X,Y,Z,omega,alpha_0,beta_0,strength_in,phi_in,theta_in,vflag);
end

WF(iiii).vals = si;
if ny == 1
    WF(iiii).vi = squeeze(sum(vi,2))';
else
    WF(iiii).vi = vi;
end

if sflag
    WF(iiii).misfit = 1;
else
    D_in = DT(iiii).vals;
    WF(iiii).misfit = 1/2.*(sum((si-D_in).^2)); 
end

end
timestep = toc(timeref);


end

