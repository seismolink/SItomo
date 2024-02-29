function [mf,G]=SIsensitivityBFGS3D(m_in,funinp,fig)
% Copyright 2024 F.Link and M.D.Long 

data = funinp.data;
indoi = funinp.indoi;
w = funinp.weights;
sflag = funinp.sflag;
vflag = funinp.vflag;
dec = funinp.dec;
thetaflag = funinp.thetaflag;
zlim = funinp.zlim;
nx = funinp.nx;
ny = funinp.ny;
nz = funinp.nz;
ix = funinp.ix;
iy = funinp.iy;
iz = funinp.iz;
edgeadd = funinp.edge;
sets = funinp.sets;
smoothpar = funinp.smooth;
initflag = funinp.initflag;
if initflag
indx = funinp.indx;
end


AZIMUTH=data.baz(indoi);
per=data.per(indoi);
par=data.p(indoi);

X_ST=data.x.*1000;
Y_ST=data.y.*1000;

X=abs(max(X_ST)-min(X_ST)+edgeadd);
dx=X/nx;
xm=linspace(min(X_ST)-(edgeadd/2),max(X_ST)+(edgeadd/2),nx);
ym=linspace(min(Y_ST)-(edgeadd/2),max(Y_ST)+(edgeadd/2),ny);
zm=linspace(zlim(1),zlim(2),nz);
X_ST = X_ST(indoi);
Y_ST = Y_ST(indoi);

if ix~=nx || iy~=ny || iz~=nz
    x=linspace(min(data.x)-(edgeadd/1000/2),max(data.x)+(edgeadd/1000/2),nx);
    y=linspace(min(data.y)-(edgeadd/1000/2),max(data.y)+(edgeadd/1000/2),ny);
    z=linspace(zlim(1),zlim(2),nz);
    xn=linspace(min(x),max(x),ix+1);
    yn=linspace(min(y),max(y),iy+1);
    zn=linspace(min(z),max(z),iz+1);
    dxn = xn(2)-xn(1);
    dyn = yn(2)-yn(1);
    dzn = zn(2)-zn(1);
    xn(end) = [];
    yn(end) = [];
    zn(end) = [];
    xn = xn+dxn/2;
    yn = yn+dyn/2;
    zn = zn+dzn/2;
    [Y,X,Z] = meshgrid(y,x,z);
    [Yn,Xn,Zn] = meshgrid(yn,xn,zn);
    k = dsearchn([Xn(:) Yn(:) Zn(:)],[X(:) Y(:) Z(:)]);
    ku = unique(k);
    if initflag
        a0 = zeros(1,ix*iy*iz);
        b0 = a0;
        c0 = a0;
        a0(indx) = m_in(1:sum(indx));
        b0(indx) = m_in(sum(indx)+1:2*sum(indx));
        c0(indx) = m_in(2*sum(indx)+1:3*sum(indx));
        strength_in = reshape(a0(k)',nx,ny,nz);
        phi_in = reshape(b0(k)',nx,ny,nz);
        theta_in = reshape(c0(k)',nx,ny,nz);
    else
        strength_in = reshape(m_in(k)',nx,ny,nz);
        phi_in = reshape(m_in(k+max(k(:)))',nx,ny,nz);
        theta_in = reshape(m_in(k+2.*max(k(:)))',nx,ny,nz);
    end
else
    if initflag
        a0 = zeros(1,nx*nz);
        b0 = a0;
        c0 = a0;
        a0(indx) = m_in(1:sum(indx));
        b0(indx) = m_in(sum(indx)+1:2*sum(indx));
        c0(indx) = m_in(2*sum(indx)+1:3*sum(indx));
        strength_in = reshape(a0,nx,ny,nz);
        phi_in = reshape(b0,nx,ny,nz);
        theta_in = reshape(c0,nx,ny,nz);
    else
        strength_in = reshape(m_in(1:nx*ny*nz),nx,ny,nz);
        phi_in = reshape(m_in(nx*ny*nz+1:2*nx*ny*nz),nx,ny,nz);
        theta_in = reshape(m_in(2*nx*ny*nz+1:3*nx*ny*nz),nx,ny,nz);
    end
end

if ~thetaflag
    theta_in = zeros(size(theta_in));
end

ioi = strength_in<0;
strength_in(ioi) = abs(strength_in(ioi));
phi_in(ioi) = mod(phi_in(ioi)+pi/2,pi);


% calculate smoothing
sig1 = strength_in.*cos(phi_in).*cos(theta_in);
sig2 = strength_in.*sin(phi_in).*cos(theta_in);
sig3 = strength_in.*sin(theta_in);
L1 = zeros(size(sig1));
L1(2:end-1,:,2:end-1) = sig1(1:end-2,:,2:end-1)+sig1(3:end,:,2:end-1)+sig1(2:end-1,:,1:end-2)+sig1(2:end-1,:,3:end)-4.*sig1(2:end-1,:,2:end-1);
L2 = zeros(size(sig2));
L2(2:end-1,:,2:end-1) = sig2(1:end-2,:,2:end-1)+sig2(3:end,:,2:end-1)+sig2(2:end-1,:,1:end-2)+sig2(2:end-1,:,3:end)-4.*sig2(2:end-1,:,2:end-1);
L3 = zeros(size(sig2));
L3(2:end-1,:,2:end-1) = sig3(1:end-2,:,2:end-1)+sig3(3:end,:,2:end-1)+sig3(2:end-1,:,1:end-2)+sig3(2:end-1,:,3:end)-4.*sig3(2:end-1,:,2:end-1);
dxsig1 = cos(phi_in).*cos(theta_in);
dxsig2 = sin(phi_in).*cos(theta_in);
dxsig3 = sin(theta_in);
dpsig1 = -sin(phi_in).*cos(theta_in).*(strength_in~=0).*sign(strength_in);
dpsig2 = cos(phi_in).*cos(theta_in).*(strength_in~=0).*sign(strength_in);
dpsig3 = zeros(size(strength_in));
dtsig1 = -cos(phi_in).*sin(theta_in).*(strength_in~=0).*sign(strength_in);
dtsig2 = -sin(phi_in).*sin(theta_in).*(strength_in~=0).*sign(strength_in);
dtsig3 = cos(theta_in).*(strength_in~=0).*sign(strength_in);
Lx = -(L1.*dxsig1+L2.*dxsig2+L3.*dxsig3);
Lp = -(L1.*dpsig1+L2.*dpsig2+L3.*dpsig3);
Lt = -(L1.*dtsig1+L2.*dtsig2+L3.*dtsig3);

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
mf = 0;

M1=max(size(xm(:)));
M2=max(size(ym(:)));
M3=max(size(zm(:)));

vpe=coeff_p(1)+coeff_p(2).*zm; %smoothed out linear P-wave speed
vse=coeff_s(1)+coeff_s(2).*zm; %smoothed out linear S-wave speed and so on

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
[Y,X,Z] = meshgrid(ym,xm,zm);

% initialize variable
wnorm = 0;
for i = 1:Nst
    WF(i).vals = 0;
    WF(i).misfit = 0;
    DT(i).vals = data.si(indoi(i));
    DT(i).err = data.err(indoi(i));
    DT(i).w = max(0,1-DT(i).err./1);
    wnorm = wnorm+DT(i).w;
end

timeref = tic;
if isempty(gcp)
    parpool(feature('NumCores'));
end
% parfor iiii = 1:Nst
for iiii = 1:Nst

% station location
x_1=X_ST(iiii);
y_1=Y_ST(iiii);
z_1=0;

[si,dSKx,dSKphi,dSKtheta] = ForwardSIsensitivity(per(iiii),par(iiii),AZIMUTH(iiii),x_1,y_1,z_1,X,Y,Z,omega,alpha_0,beta_0,strength_in,phi_in,theta_in,dec,vflag);

WF(iiii).vals = si;

if sflag
    WF(iiii).misfit = 1;
else
    D_in = DT(iiii).vals;
    WF(iiii).misfit = si-D_in;
end
mf = mf+1/2.*(WF(iiii).misfit.^2).*DT(iiii).w;

temp1a = real(dSKx);
temp2a = real(dSKphi);
temp3a = real(dSKtheta);
temp1 = temp1a;
temp2 = temp2a;
temp3 = temp3a;
if~thetaflag
    temp3 = zeros(size(temp3));
end
EKx = EKx+temp1.*w(iiii).*WF(iiii).misfit.*DT(iiii).w;
EKphi = EKphi+temp2.*w(iiii).*WF(iiii).misfit.*DT(iiii).w;
EKtheta = EKtheta+temp3.*w(iiii).*WF(iiii).misfit.*DT(iiii).w;
end
EKxa = transform_kernel3D(EKx,sets);
EKphia = transform_kernel3D(EKphi,sets);
EKthetaa = transform_kernel3D(EKtheta,sets);

Rx = -(smoothpar).*max(abs(EKxa(:)))./max(abs(Lx(:))).*Lx;
Rx(isnan(Rx)) = 0; Rx(isinf(Rx)) = 0;
Rp = -(smoothpar).*max(abs(EKphia(:)))./max(abs(Lp(:))).*Lp;
Rp(isnan(Rp)) = 0; Rp(isinf(Rp)) = 0;
Rt = -(smoothpar).*max(abs(EKthetaa(:)))./max(abs(Lt(:))).*Lt;
Rt(isnan(Rt)) = 0; Rt(isinf(Rt)) = 0;
EKx = EKxa+Rx;
EKphi = EKphia+Rp;
EKtheta = EKthetaa+Rt;

if ix~=nx || iy~=ny || iz~=nz
    EKx = squeeze(EKx);
    EKphi = squeeze(EKphi);
    EKtheta = squeeze(EKtheta);
    for i = 1:length(ku)
        EKxn(i) = mean(EKx(k==ku(i)));
        EKphin(i) = mean(EKphi(k==ku(i)));
        EKthetan(i) = mean(EKtheta(k==ku(i)));
    end
    EKx = reshape(EKxn,ix,iy,iz);
    EKphi = reshape(EKphin,ix,iy,iz);
    EKtheta = reshape(EKthetan,ix,iy,iz);
end
if initflag
    a0 = -EKx(:);
    b0 = -EKphi(:);
    c0 = EKtheta(:);
    G = [a0(indx)' b0(indx)' c0(indx)'];
else
    G = [-EKx(:)' -EKphi(:)' EKtheta(:)']; 
end
G = G./max(abs(G(:)));

mf2 = 1/2.*sum((L1(:)).^2+(L2(:)).^2+(L3(:)).^2)*length(WF)/4;
mf = mf+mf2.*smoothpar.^2;

    xx = unique(data.x(indoi));
    yy = unique(data.y(indoi));
    [XX,YY] = meshgrid(xx,yy);
    XX = XX';
    YY = YY';
    sioi = zeros(length(xx),length(yy));
    coi = zeros(size(sioi));
    for i = 1:length(indoi)
        dd1 = sqrt((xx-data.x(indoi(i))).^2);
        dd2 = sqrt((yy-data.y(indoi(i))).^2);
        [~,I] = min(dd1);
        [~,J] = min(dd2);
        sioi(I,J) = sioi(I,J)+1/2.*WF(i).misfit.^2;
        coi(I,J) = coi(I,J)+1;
    end
    sioi = sioi./coi;
    sioi(coi==0)=0;
    sipl = sioi(coi>0);
    cpl = coi(coi>0);
    cpl = cpl-min(cpl);
    cpl = cpl./(0.3*max(cpl))+1;
    ypl = YY(coi>0);
    xpl = XX(coi>0);
    silim = [0 max(abs(sioi(:)))];
    sipl0 = sipl;
    xpl0 = xpl-(xx(2)-xx(1))/10;
    ypl0 = ypl-(yy(2)-yy(1))/10;
    
    if funinp.showflag
        figure(fig);
    else
        set(0, 'CurrentFigure', fig);
    end
    subplot(4,1,1)
    scatter(xpl0,ypl0,cpl*30,sipl0,'filled');
    hold on;
    scatter(xpl,ypl,cpl*30,sipl,'filled');
    axis([min(xpl)-5*(xx(2)-xx(1))/10 max(xpl)+5*(xx(2)-xx(1))/10 min(ypl)-5*(yy(2)-yy(1))/10 max(ypl)+5*(yy(2)-yy(1))/10 ])
    ylabel('y in [km]','FontSize',9)
    caxis(silim)
    cax = colorbar;
    ylabel(cax,'Misfit','FontSize',9)
    hold off
    subplot(4,1,2)
    imagesc(x,z,mod(squeeze(phi_in(:,round(min(iy,ny)/2),:).*180./pi),180)')
    axis([min(xpl)-5*(xx(2)-xx(1))/10 max(xpl)+5*(xx(2)-xx(1))/10 min(z) max(z)])
    ylabel('z in [km]','FontSize',9)
    caxis([0 180]);
    cax1 = colorbar;
    ylabel(cax1,'Phi in [deg]','FontSize',9)
    subplot(4,1,3)
    imagesc(x,z,squeeze(strength_in(:,round(min(iy,ny)/2),:)*100)')
    axis([min(xpl)-5*(xx(2)-xx(1))/10 max(xpl)+5*(xx(2)-xx(1))/10 min(z) max(z)])
    ylabel('z in [km]','FontSize',9)
    caxis([0 100]);
    caxis([-1 25])
    cax2 = colorbar;
    ylabel(cax2,'Ratio in [pct]','FontSize',9)
    subplot(4,1,4)
    imagesc(x,z,squeeze(theta_in(:,round(min(iy,ny)/2),:)*180/pi)')
    axis([min(xpl)-5*(xx(2)-xx(1))/10 max(xpl)+5*(xx(2)-xx(1))/10 min(z) max(z)])
    ylabel('z in [km]','FontSize',9)
    xlabel('x in [km]','FontSize',9)
    caxis([-60 60])
    cax3 = colorbar;
    ylabel(cax3,'Theta in [deg]','FontSize',9)
    drawnow
if ~funinp.showflag
    fid = fopen(funinp.filename,'at');
    fprintf(fid,'Current misfit = %s.\n',mf);
    fclose(fid);
else
    disp(['Current misfit = ' num2str(mf)])
end

end

