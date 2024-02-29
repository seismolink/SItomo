function [model,result] = DoInvert3Dmcmc(data,thetaflag,vflag,zlim,nx,ny,nz,edgeadd,xstepini,phistepini,thetastepini,N1,N2,n1,n2,nevmax,initflag,strength_in,phi_in,theta_in,smoothpar,showflag,printflag,graphicsfolder,filename,addname)
% Copyright 2024 F.Link and M.D.Long 

[workdir,~,~] = fileparts(graphicsfolder(1:end-1));
if ~showflag
    fid = fopen([workdir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Initialize results variables.\n');
    fclose(fid);
else
    disp('Initialize results variables.')
end

if ~showflag
    fid = fopen([workdir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Start loop for statistics.\n');
    fclose(fid);
else
    disp('Start loop for statistics.')
end

A_0 = zeros(nx,ny,nz);
B_0 = zeros(nx,ny,nz);
C_0 = zeros(nx,ny,nz);

xstep = xstepini;
phistep = phistepini;
thetastep = thetastepini;

% define positionvectors
X_ST=data.x;
Y_ST=data.y;
x=linspace(min(X_ST)-(edgeadd/1000/2),max(X_ST)+(edgeadd/1000/2),nx);
y=linspace(min(Y_ST)-(edgeadd/1000/2),max(Y_ST)+(edgeadd/1000/2),ny);
z=linspace(zlim(1),zlim(2),nz);
[Y,X,Z] = meshgrid(y,x,z);
    
% initial Voronoi cells
noc = N1;
ifpmts = 1:length(X(:));
xini = 0.0.*ones(noc,1);
phiini = zeros(noc,1);
thetaini = zeros(noc,1);
vn = randperm(length(X(:)),noc);
rjv = zeros(size(vn));
xv = X(vn);
yv = Y(vn);
zv = Z(vn);
ifpmts(vn) = [];
Idx = knnsearch([xv(:) yv(:) zv(:)],[X(:) Y(:) Z(:)]);
for i = 1:noc
    A_0(Idx==i) = xini(i);
    B_0(Idx==i) = phiini(i);
    C_0(Idx==i) = thetaini(i);
end

% stepsize for nuclei movement (standarddeviation sigma = 2*dx
xu = unique(X(:));
ddx = xu(2)-xu(1);
sigma = 3*ddx;

% calculate initial splitting intensities
if ~showflag
    fid = fopen([workdir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Calculate reference SI.\n');
    fclose(fid);
else
    disp('Calculate reference SI.')
end
[WFref,~]=SIforward(data,1:length(data.si),0,vflag,A_0,B_0,C_0,zlim.*1000,nx,ny,nz,edgeadd,1);
for i = 1:length(WFref)
    siref(i) = WFref(i).vals;
end

% Prepare shape of likelyhood function
roi = 0;
sig = 0;
while roi<0.01%eps
    sig = sig+0.01;
    L1 = likelyhood_smooth(data,1:length(data.si),data.si,sig,A_0,B_0,C_0,smoothpar);
    L2 = likelyhood_smooth(data,1:length(data.si),siref,sig,A_0,B_0,C_0,smoothpar);
    roi = L2./L1;
end

% calculate reference likelyhood
Lold = likelyhood_smooth(data,1:length(data.si),siref,sig,A_0,B_0,C_0,smoothpar);
noi = min(nevmax,round(length(data.x)*2/4));

Nmc = N2;
nfail = 0;
nfail2 = 0;
if printflag
    if showflag
        fig = figure('Position',[300 50 600 580]);
    else
        fig = figure('Position',[300 50 600 580],'visible','off');
    end
    splt1 = subplot(3,2,[1 3 5]);
    hold on
    splt2 = subplot(3,2,2);
    splt3 = subplot(3,2,4);
    splt4 = subplot(3,2,6);
end


indoi = randperm(length(data.x),noi);
for nmc = 1:Nmc
    if nmc < 5
        jj = noc;
    else
        jj = fix(noc./20);
    end
% Recalculate reference
    if ~mod(nmc,n1)
        jj = noc;
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Re-calculate reference SI.\n');
            fclose(fid);
        else
            disp('Re-calculate reference SI.')
        end
        indoio = indoi;
        indoi = randperm(length(data.x),noi);
        [WFrefn,~]=SIforward(data,1:length(data.si),0,vflag,A_0,B_0,C_0,zlim.*1000,nx,ny,nz,edgeadd,1);
        dref = 0;
        for i = 1:length(WFrefn)
            if ~isempty(find(indoio==i,1))
                dref = dref+sqrt((WFrefn(i).vals-siref(i)).^2);
            end
            siref(i) = WFrefn(i).vals;
        end
        WFref = WFrefn;
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Corrected for SIdiff = %f.\n',dref);
            fclose(fid);
        else
            fprintf('Corrected for SIdiff = %f.\n',dref);
        end
        Lold = likelyhood_smooth(data,indoi,siref(indoi),sig,A_0,B_0,C_0,smoothpar);
    end
    xv = X(vn);
    yv = Y(vn);
    zv = Z(vn);
    Idx = knnsearch([xv(:) yv(:) zv(:)],[X(:) Y(:) Z(:)]);
        
    % draw subset of data
    
    if ~mod(nmc/2,2)
        thetaflag2 = 1;
    else
        thetaflag2 = 0;
    end
    if ~mod(nmc,2)
        % even step >> Perturbation of the model
        % --------------------------------------
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'%f of %f - Even step >> Perturbation of the model.\n', nmc,Nmc);
            fclose(fid);
        else
            fprintf('%f of %f - Even step >> Perturbation of the model.\n', nmc,Nmc);
        end
        for ii = 1:jj
        % select a voronoi cell
        ioi = randi(noc,1);
        nt = 1;
        nm = 30;
        while rjv(ioi)>min(rjv)&&nt<nm
            ioi = randi(noc,1);
            nt = nt+1;
        end
        if sum(Idx==ioi) == 0
            if ~showflag
                fid = fopen([workdir '/results/log' addname '.txt'],'at');
                fprintf(fid,'ERROR: Problem with voronoi-cell selection.\n');
                fclose(fid);
                return
            else
                disp('ERROR: Problem with voronoi-cell selection.\n');
                return
            end
        end
        % calculate gradient for phi
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Calculate reference Gradient for phi.\n');
            fclose(fid);
        else
            disp('Calculate reference Gradient for phi.')
        end
        [EKphi]=phisensitivity(data,WFref(indoi),indoi,vflag,A_0,B_0,C_0,Idx==ioi,zlim.*1000,nx,ny,nz,edgeadd,1);
        
        sig1 = A_0.*cos(B_0).*cos(C_0);
        sig2 = A_0.*sin(B_0).*cos(C_0);
        sig3 = A_0.*sin(C_0);
        L1 = zeros(size(sig1));
        L1(2:end-1,:,2:end-1) = sig1(1:end-2,:,2:end-1)+sig1(3:end,:,2:end-1)+sig1(2:end-1,:,1:end-2)+sig1(2:end-1,:,3:end)-4.*sig1(2:end-1,:,2:end-1);
        L2 = zeros(size(sig2));
        L2(2:end-1,:,2:end-1) = sig2(1:end-2,:,2:end-1)+sig2(3:end,:,2:end-1)+sig2(2:end-1,:,1:end-2)+sig2(2:end-1,:,3:end)-4.*sig2(2:end-1,:,2:end-1);
        L3 = zeros(size(sig2));
        L3(2:end-1,:,2:end-1) = sig3(1:end-2,:,2:end-1)+sig3(3:end,:,2:end-1)+sig3(2:end-1,:,1:end-2)+sig3(2:end-1,:,3:end)-4.*sig3(2:end-1,:,2:end-1);
        dpsig1 = -sin(B_0).*cos(C_0).*(A_0~=0).*sign(A_0);
        dpsig2 = cos(B_0).*cos(C_0).*(A_0~=0).*sign(A_0);
        dpsig3 = zeros(size(A_0));
        Lp = -(L1.*dpsig1+L2.*dpsig2+L3.*dpsig3);
        Rp = -(smoothpar).*max(abs(EKphi(:)))./max(abs(Lp(:))).*Lp;
        Rp(isnan(Rp)) = 0; Rp(isinf(Rp)) = 0;
        EKphi = EKphi+mean(Rp(Idx==ioi));
        
        id = rand(10,1);
        % random step in phi
        dB = phistep.*id(randi(10,1));
        if EKphi<0
            dB = -dB;
        end
        B0 = squeeze(B_0);
        B0(Idx==ioi) = mod(B0(Idx==ioi)-dB,pi);
        % calculate gradient for x
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Calculate reference Gradient for X.\n');
            fclose(fid);
        else
            disp('Calculate reference Gradient for X.')
        end
        [EKx]=xsensitivity(data,WFref(indoi),indoi,vflag,A_0,B_0,C_0,Idx==ioi,zlim.*1000,nx,ny,nz,edgeadd,1);
        
        sig1 = A_0.*cos(B_0).*cos(C_0);
        sig2 = A_0.*sin(B_0).*cos(C_0);
        sig3 = A_0.*sin(C_0);
        L1 = zeros(size(sig1));
        L1(2:end-1,:,2:end-1) = sig1(1:end-2,:,2:end-1)+sig1(3:end,:,2:end-1)+sig1(2:end-1,:,1:end-2)+sig1(2:end-1,:,3:end)-4.*sig1(2:end-1,:,2:end-1);
        L2 = zeros(size(sig2));
        L2(2:end-1,:,2:end-1) = sig2(1:end-2,:,2:end-1)+sig2(3:end,:,2:end-1)+sig2(2:end-1,:,1:end-2)+sig2(2:end-1,:,3:end)-4.*sig2(2:end-1,:,2:end-1);
        L3 = zeros(size(sig2));
        L3(2:end-1,:,2:end-1) = sig3(1:end-2,:,2:end-1)+sig3(3:end,:,2:end-1)+sig3(2:end-1,:,1:end-2)+sig3(2:end-1,:,3:end)-4.*sig3(2:end-1,:,2:end-1);
        dxsig1 = cos(B_0).*cos(C_0);
        dxsig2 = sin(B_0).*cos(C_0);
        dxsig3 = sin(C_0);
        Lx = -(L1.*dxsig1+L2.*dxsig2+L3.*dxsig3);
        Rx = -(smoothpar).*max(abs(EKx(:)))./max(abs(Lx(:))).*Lx;
        Rx(isnan(Rx)) = 0; Rx(isinf(Rx)) = 0;
        EKx = EKx+mean(Rx(Idx==ioi));
        % random step in x
        dA = xstep.*id(randi(10,1));
        if EKx<0
            dA = -dA;
        end
        A0 = squeeze(A_0);
        A0(Idx==ioi) = A0(Idx==ioi)'-dA;
        A0(A0>1) = 1;
        A0(A0<0) = 0;
        % random step in theta
        if thetaflag && thetaflag2
        % calculate gradient for theta
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Calculate reference Gradient for theta.\n');
            fclose(fid);
        else
            disp('Calculate reference Gradient for theta.')
        end
        [EKtheta]=thetasensitivity(data,WFref(indoi),indoi,vflag,A_0,B_0,C_0,Idx==ioi,zlim.*1000,nx,ny,nz,edgeadd,1);
        
        sig1 = A_0.*cos(B_0).*cos(C_0);
        sig2 = A_0.*sin(B_0).*cos(C_0);
        sig3 = A_0.*sin(C_0);
        L1 = zeros(size(sig1));
        L1(2:end-1,:,2:end-1) = sig1(1:end-2,:,2:end-1)+sig1(3:end,:,2:end-1)+sig1(2:end-1,:,1:end-2)+sig1(2:end-1,:,3:end)-4.*sig1(2:end-1,:,2:end-1);
        L2 = zeros(size(sig2));
        L2(2:end-1,:,2:end-1) = sig2(1:end-2,:,2:end-1)+sig2(3:end,:,2:end-1)+sig2(2:end-1,:,1:end-2)+sig2(2:end-1,:,3:end)-4.*sig2(2:end-1,:,2:end-1);
        L3 = zeros(size(sig2));
        L3(2:end-1,:,2:end-1) = sig3(1:end-2,:,2:end-1)+sig3(3:end,:,2:end-1)+sig3(2:end-1,:,1:end-2)+sig3(2:end-1,:,3:end)-4.*sig3(2:end-1,:,2:end-1);
        dtsig1 = -cos(B_0).*sin(C_0).*(A_0~=0).*sign(A_0);
        dtsig2 = -sin(B_0).*sin(C_0).*(A_0~=0).*sign(A_0);
        dtsig3 = cos(C_0).*(A_0~=0).*sign(A_0);
        Lt = -(L1.*dtsig1+L2.*dtsig2+L3.*dtsig3);
        Rt = -(smoothpar).*max(abs(EKtheta(:)))./max(abs(Lt(:))).*Lt;
        Rt(isnan(Rt)) = 0; Rt(isinf(Rt)) = 0;
        EKtheta = EKtheta+mean(Rt(Idx==ioi));
        
            dC = xstep.*id;
            if EKtheta>0
                dC = -dC;
            end
            C0 = squeeze(C_0);
            C0(Idx==ioi) = C_0(Idx==ioi)-dC(randi(10,1));
            C0(C0>pi/3) = pi/3;
            C0(C0<-pi/3) = -pi/3;
        else
            dC = 0;
            C0 = squeeze(C_0);
        end
        A0 = reshape(A0,nx,ny,nz);
        B0 = reshape(B0,nx,ny,nz);
        C0 = reshape(C0,nx,ny,nz);
        % calculate change in splitting intensity
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Calculate change in splitting intensity.\n');
            fclose(fid);
        else
            disp('Calculate change in splitting intensity.')
        end
        [dWF]=dsiforward(data,indoi,vflag,A0,B0,C0,Idx==ioi,zlim.*1000,nx,ny,nz,edgeadd,1);
        for i = 1:length(dWF)
            si(i) = siref(indoi(i))+sum(dWF(i).vi'-WFref(indoi(i)).vi(Idx==ioi));
        end
        L = likelyhood_smooth(data,indoi,si,sig,A0,B0,C0,smoothpar);
        prorat = 1;
        % acceptance probability
        alpha = L/Lold*prorat;        
        crit = 1; % try to enforce improvement
        nc = 3;
        nt = 1;
        while alpha<crit && nt < nc
            dA = dA/2;
            dB = dB/2;
            dC = dC/2;
            A0 = squeeze(A_0);
            A0(Idx==ioi) = A0(Idx==ioi)'-dA;
            A0(A0>1) = 1;
            A0(A0<0) = 0;
            B0 = squeeze(B_0);
            B0(Idx==ioi) = mod(B0(Idx==ioi)-dB,pi);
            C0 = squeeze(C_0);
            C0(Idx==ioi) = C_0(Idx==ioi)-dC;
            C0(C0>pi/3) = pi/3;
            C0(C0<-pi/3) = -pi/3;
            A0 = reshape(A0,nx,ny,nz);
            B0 = reshape(B0,nx,ny,nz);
            C0 = reshape(C0,nx,ny,nz);
            
            [dWF]=dsiforward(data,indoi,vflag,A0,B0,C0,Idx==ioi,zlim.*1000,nx,ny,nz,edgeadd,1);
            for i = 1:length(dWF)
                si(i) = siref(indoi(i))+sum(dWF(i).vi'-WFref(indoi(i)).vi(Idx==ioi));
            end
            L = likelyhood_smooth(data,indoi,si,sig,A0,B0,C0,smoothpar);
            % acceptance probability
            alpha = L/Lold*prorat;  
            nt = nt+1;
        end
        crit = 0.99995+rand(1)*0.00005;
        if alpha>=crit
            if ~showflag
                fid = fopen([workdir '/results/log' addname '.txt'],'at');
                fprintf(fid,'Accept current step.\n');
                fprintf(fid,'alpha = %f | crit = %f | LH = %f\n',alpha,crit,L);
                fclose(fid);
            else
                disp('Accept current step.')
                fprintf('alpha = %f | crit = %f | LH = %f\n',alpha,crit,L);
            end
            A_0 = A0;
            B_0 = B0;
            C_0 = C0;
            for i = 1:length(si)
                WFref(indoi(i)).vi(Idx==ioi) = dWF(i).vi;%(Idx==ioi);
                WFref(indoi(i)).vals = si(i);
            end
            MCMC(nmc).L = L;
            MCMC(nmc).x = A0(:);
            MCMC(nmc).phi = B0(:);
            MCMC(nmc).theta = C0(:);
            A0 = squeeze(A_0);
            B0 = squeeze(B_0);
            C0 = squeeze(C_0);
            MCMC(nmc).a0 = A0(newvn);
            MCMC(nmc).b0 = B0(newvn);
            MCMC(nmc).c0 = C0(newvn);
            MCMC(nmc).vn = vn;
            MCMC(nmc).xv = X(vn);
            MCMC(nmc).yv = Y(vn);
            MCMC(nmc).zv = Z(vn);
            xv = X(vn);
            yv = Y(vn);
            zv = Z(vn);
            Idx = knnsearch([xv(:) yv(:) zv(:)],[X(:) Y(:) Z(:)]);
            Lold = L;
            nfail = 0;
            siref(indoi) = si;
        else
            if ~showflag
                fid = fopen([workdir '/results/log' addname '.txt'],'at');
                fprintf(fid,'Reject current step.\n');
                fprintf(fid,'alpha = %f | crit = %f | LH = %f\n',alpha,crit,Lold);
                fclose(fid);
            else
                disp('Reject current step.')
                fprintf('alpha = %f | crit = %f | LH = %f\n',alpha,crit,Lold);
            end
            MCMC(nmc).L = Lold;
            MCMC(nmc).x = A_0(:);
            MCMC(nmc).phi = B_0(:);
            MCMC(nmc).theta = C_0(:);
            A0 = squeeze(A_0);
            B0 = squeeze(B_0);
            C0 = squeeze(C_0);
            MCMC(nmc).a0 = A0(newvn);
            MCMC(nmc).b0 = B0(newvn);
            MCMC(nmc).c0 = C0(newvn);
            MCMC(nmc).vn = vn;
            MCMC(nmc).xv = X(vn);
            MCMC(nmc).yv = Y(vn);
            MCMC(nmc).zv = Z(vn);
            xv = X(vn);
            yv = Y(vn);
            zv = Z(vn);
            rjv(ioi) = rjv(ioi)+1;
            Idx = knnsearch([xv(:) yv(:) zv(:)],[X(:) Y(:) Z(:)]);
            nfail = nfail+1;
        end
        end
    else
        % odd step >> perturbation of cellular parametrization (birth,
        % death, move)
        % --------------------------------------
        bdmflag = randi(3,1);
        if bdmflag == 1
            if ~showflag
                fid = fopen([workdir '/results/log' addname '.txt'],'at');
                fprintf(fid,'%f of %f - Odd step >> Birth of a new cell.\n', nmc,Nmc);
                fclose(fid);
            else
                fprintf('%f of %f - Odd step >> Birth of a new cell.\n', nmc,Nmc);
            end
            % birth of a new voronoi cell
            temp = ifpmts(randi(length(ifpmts),1));
            noc = noc+1;
            newvn = [vn temp];
            rjvn = [rjv 0];
            xv = X(newvn);
            yv = Y(newvn);
            zv = Z(newvn);
            Idxo = Idx;
            Idx = knnsearch([xv(:) yv(:) zv(:)],[X(:) Y(:) Z(:)]);
            id2 = abs(Idxo-Idx)>0;
            ioi = length(newvn);
            A0 = squeeze(A_0);
            B0 = squeeze(B_0);
            C0 = squeeze(C_0);
            a0 = A0(newvn);
            b0 = B0(newvn);
            c0 = C0(newvn);
            for i = 1:length(newvn)
                A0(Idx==i)=a0(i);
                B0(Idx==i)=b0(i);
                C0(Idx==i)=c0(i);
            end
            % calculate perturbation for new model
            
            % calculate gradient for phi
            [EKphi]=phisensitivity(data,WFref(indoi),indoi,vflag,A0,B0,C0,Idx==ioi,zlim.*1000,nx,ny,nz,edgeadd,1);
            sig1 = A0.*cos(B0).*cos(C0);
            sig2 = A0.*sin(B0).*cos(C0);
            sig3 = A0.*sin(C0);
            L1 = zeros(size(sig1));
            L1(2:end-1,:,2:end-1) = sig1(1:end-2,:,2:end-1)+sig1(3:end,:,2:end-1)+sig1(2:end-1,:,1:end-2)+sig1(2:end-1,:,3:end)-4.*sig1(2:end-1,:,2:end-1);
            L2 = zeros(size(sig2));
            L2(2:end-1,:,2:end-1) = sig2(1:end-2,:,2:end-1)+sig2(3:end,:,2:end-1)+sig2(2:end-1,:,1:end-2)+sig2(2:end-1,:,3:end)-4.*sig2(2:end-1,:,2:end-1);
            L3 = zeros(size(sig2));
            L3(2:end-1,:,2:end-1) = sig3(1:end-2,:,2:end-1)+sig3(3:end,:,2:end-1)+sig3(2:end-1,:,1:end-2)+sig3(2:end-1,:,3:end)-4.*sig3(2:end-1,:,2:end-1);
            dpsig1 = -sin(B0).*cos(C0).*(A0~=0).*sign(A0);
            dpsig2 = cos(B0).*cos(C0).*(A0~=0).*sign(A0);
            dpsig3 = zeros(size(A0));
            Lp = -(L1.*dpsig1+L2.*dpsig2+L3.*dpsig3);
            Rp = -(smoothpar).*max(abs(EKphi(:)))./max(abs(Lp(:))).*Lp;
            Rp(isnan(Rp)) = 0; Rp(isinf(Rp)) = 0;
            EKphi = EKphi+mean(Rp(Idx==ioi));
            % calculate gradient for x
            [EKx]=xsensitivity(data,WFref(indoi),indoi,vflag,A0,B0,C0,Idx==ioi,zlim.*1000,nx,ny,nz,edgeadd,1);
            sig1 = A0.*cos(B0).*cos(C0);
            sig2 = A0.*sin(B0).*cos(C0);
            sig3 = A0.*sin(C0);
            L1 = zeros(size(sig1));
            L1(2:end-1,:,2:end-1) = sig1(1:end-2,:,2:end-1)+sig1(3:end,:,2:end-1)+sig1(2:end-1,:,1:end-2)+sig1(2:end-1,:,3:end)-4.*sig1(2:end-1,:,2:end-1);
            L2 = zeros(size(sig2));
            L2(2:end-1,:,2:end-1) = sig2(1:end-2,:,2:end-1)+sig2(3:end,:,2:end-1)+sig2(2:end-1,:,1:end-2)+sig2(2:end-1,:,3:end)-4.*sig2(2:end-1,:,2:end-1);
            L3 = zeros(size(sig2));
            L3(2:end-1,:,2:end-1) = sig3(1:end-2,:,2:end-1)+sig3(3:end,:,2:end-1)+sig3(2:end-1,:,1:end-2)+sig3(2:end-1,:,3:end)-4.*sig3(2:end-1,:,2:end-1);
            dxsig1 = cos(B0).*cos(C0);
            dxsig2 = sin(B0).*cos(C0);
            dxsig3 = sin(C0);
            Lx = -(L1.*dxsig1+L2.*dxsig2+L3.*dxsig3);
            Rx = -(smoothpar).*max(abs(EKx(:)))./max(abs(Lx(:))).*Lx;
            Rx(isnan(Rx)) = 0; Rx(isinf(Rx)) = 0;
            EKx = EKx+mean(Rx(Idx==ioi));
            id = rand(10,1);
            % random step in phi
            dB = phistep.*id(randi(10,1));
            if EKphi<0
                dB = -dB;
            end
            B0(Idx==ioi) = mod(B0(Idx==ioi)-dB,pi);
            % random step in x
            dA = xstep.*id(randi(10,1));
            if EKx<0
                dA = -dA;
            end
            A0(Idx==ioi) = A0(Idx==ioi)-dA;
            A0(A0>1) = 1;
            A0(A0<0) = 0;
            % random step in theta
            if thetaflag && thetaflag2
                % calculate gradient for theta
                [EKtheta]=thetasensitivity(data,WFref(indoi),indoi,vflag,A0,B0,C0,Idx==ioi,zlim.*1000,nx,ny,nz,edgeadd,1);
                sig1 = A0.*cos(B0).*cos(C0);
                sig2 = A0.*sin(B0).*cos(C0);
                sig3 = A0.*sin(C0);
                L1 = zeros(size(sig1));
                L1(2:end-1,:,2:end-1) = sig1(1:end-2,:,2:end-1)+sig1(3:end,:,2:end-1)+sig1(2:end-1,:,1:end-2)+sig1(2:end-1,:,3:end)-4.*sig1(2:end-1,:,2:end-1);
                L2 = zeros(size(sig2));
                L2(2:end-1,:,2:end-1) = sig2(1:end-2,:,2:end-1)+sig2(3:end,:,2:end-1)+sig2(2:end-1,:,1:end-2)+sig2(2:end-1,:,3:end)-4.*sig2(2:end-1,:,2:end-1);
                L3 = zeros(size(sig2));
                L3(2:end-1,:,2:end-1) = sig3(1:end-2,:,2:end-1)+sig3(3:end,:,2:end-1)+sig3(2:end-1,:,1:end-2)+sig3(2:end-1,:,3:end)-4.*sig3(2:end-1,:,2:end-1);
                dtsig1 = -cos(B0).*sin(C0).*(A0~=0).*sign(A0);
                dtsig2 = -sin(B0).*sin(C0).*(A0~=0).*sign(A0);
                dtsig3 = cos(C0).*(A0~=0).*sign(A0);
                Lt = -(L1.*dtsig1+L2.*dtsig2+L3.*dtsig3);
                Rt = -(smoothpar).*max(abs(EKtheta(:)))./max(abs(Lt(:))).*Lt;
                Rt(isnan(Rt)) = 0; Rt(isinf(Rt)) = 0;
                EKtheta = EKtheta+mean(Rt(Idx==ioi));
            
                dC = xstep.*id;
                if EKtheta>0
                    dC = -dC;
                end
                C0(Idx==ioi) = C0(Idx==ioi)-dC(randi(10,1));
                C0(C0>pi/3) = pi/3;
                C0(C0<-pi/3) = -pi/3;
            end
            A0 = reshape(A0,nx,ny,nz);
            B0 = reshape(B0,nx,ny,nz);
            C0 = reshape(C0,nx,ny,nz);
            % calculate change in splitting intensity
            if isempty(id2)
                si = siref(indoi);
            else
                [dWF]=dsiforward(data,indoi,vflag,A0,B0,C0,id2,zlim.*1000,nx,ny,nz,edgeadd,1);
                for i = 1:length(dWF)
                    si(i) = siref(indoi(i))+sum(dWF(i).vi'-WFref(indoi(i)).vi(id2));
                end
            end
            L = likelyhood_smooth(data,indoi,si,sig,A0,B0,C0,smoothpar);
            prorat = 1;
            ioi = id2;
        end
        if bdmflag == 2
            if length(vn) <= n2
                L = 0;
                prorat = 1;
            else
            % death of a voronoi cell
            rjvn = rjv;
            if ~showflag
                fid = fopen([workdir '/results/log' addname '.txt'],'at');
                fprintf(fid,'%f of %f - Odd step >> Death of a cell.\n', nmc,Nmc);
                fclose(fid);
            else
                fprintf('%f of %f - Odd step >> Death of a cell.\n', nmc,Nmc);
            end
            temp = vn(randi(length(vn),1));
            noc = noc-1;
            newvn = vn;
            newvn(vn==temp) = [];
            rjvn(vn==temp) = [];
            xv = X(newvn);
            yv = Y(newvn);
            zv = Z(newvn);
            Idxo = Idx;
            Idx = knnsearch([xv(:) yv(:) zv(:)],[X(:) Y(:) Z(:)]);
            Idx2 = Idx;
            Idx(Idx>=find(vn==temp,1)) = Idx(Idx>=find(vn==temp,1))+1;
            id2 = abs(Idxo-Idx)>0;
            A0 = squeeze(A_0);
            B0 = squeeze(B_0);
            C0 = squeeze(C_0);
            a0 = A0(newvn);
            b0 = B0(newvn);
            c0 = C0(newvn);
            for i = 1:length(newvn)
                A0(Idx2==i)=a0(i);
                B0(Idx2==i)=b0(i);
                C0(Idx2==i)=c0(i);
            end    
            A0 = reshape(A0,nx,ny,nz);
            B0 = reshape(B0,nx,ny,nz);
            C0 = reshape(C0,nx,ny,nz);
            % calculate si correction for new model
            [dWF]=dsiforward(data,indoi,vflag,A0,B0,C0,id2,zlim.*1000,nx,ny,nz,edgeadd,1);
            for i = 1:length(dWF)
                si(i) = siref(indoi(i))+sum(dWF(i).vi'-WFref(indoi(i)).vi(id2));
            end
            L = likelyhood_smooth(data,indoi,si,sig,A0,B0,C0,smoothpar);
            prorat = 1;
            ioi = id2;
            end
        end
        if bdmflag == 3
            rjvn = rjv;
            % move of a voronoi cell
            if ~showflag
                fid = fopen([workdir '/results/log' addname '.txt'],'at');
                fprintf(fid,'%f of %f - Odd step >> Move of a cell.\n', nmc,Nmc);
                fclose(fid);
            else
                fprintf('%f of %f - Odd step >> Move of a cell.\n', nmc,Nmc);
            end
            id = 0;
            while sum(id) == 0
                temp1 = vn(randi(length(vn),1));
                x1 = xv((vn==temp1))+randn.*sigma;
                y1 = yv((vn==temp1))+randn.*sigma;
                z1 = zv((vn==temp1))+randn.*sigma;
                dd = sqrt((X(ifpmts)-x1).^2+(Y(ifpmts)-y1).^2+(Z(ifpmts)-z1).^2);
                [~,ii] = min(dd);
                temp2 = ifpmts(ii);
                newvn = vn;
                newvn(vn==temp1) = temp2;
                ifpmtso = ifpmts;
                ifpmts = sort([ifpmts vn(vn==temp1)]);
                xv = X(newvn);
                yv = Y(newvn);
                zv = Z(newvn);
                Idxo = Idx;
                Idx = knnsearch([xv(:) yv(:) zv(:)],[X(:) Y(:) Z(:)]);
                id = abs(Idx-Idxo)>0;
                if sum(id) == 0
                    Idx = Idxo;
                    ifpmts = ifpmtso;
                end
            end
            A0 = squeeze(A_0);
            B0 = squeeze(B_0);
            C0 = squeeze(C_0);
            a0 = A0(newvn);
            b0 = B0(newvn);
            c0 = C0(newvn);
            for i = 1:length(newvn)
                A0(Idx==i)=a0(i);
                B0(Idx==i)=b0(i);
                C0(Idx==i)=c0(i);
            end
            A0 = reshape(A0,nx,ny,nz);
            B0 = reshape(B0,nx,ny,nz);
            C0 = reshape(C0,nx,ny,nz);
            % calculate si correction for new model
            [dWF]=dsiforward(data,indoi,vflag,A0,B0,C0,id,zlim.*1000,nx,ny,nz,edgeadd,1);
            for i = 1:length(dWF)
                si(i) = siref(indoi(i))+sum(dWF(i).vi'-WFref(indoi(i)).vi(id));
            end
            L = likelyhood_smooth(data,indoi,si,sig,A0,B0,C0,smoothpar);
            
            prorat = 1;
            ioi = id;
        end
        
        
        % acceptance probability
        alpha = L/Lold*prorat;
        crit = 0.99995+rand(1)*0.00005;
        if alpha>=crit
            if ~showflag
                fid = fopen([workdir '/results/log' addname '.txt'],'at');
                fprintf(fid,'Accept current step.\n');
                fprintf(fid,'alpha = %f | crit = %f | LH = %f\n',alpha,crit,L);
                fclose(fid);
            else
                disp('Accept current step.')
                fprintf('alpha = %f | crit = %f | LH = %f\n',alpha,crit,L);
            end
            A_0 = A0;
            B_0 = B0;
            C_0 = C0;
            vn = newvn;
            rjv = rjvn;
            noc = length(vn);
            for i = 1:length(si)
                WFref(indoi(i)).vi(ioi) = dWF(i).vi;
                WFref(indoi(i)).vals = si(i);
            end
            MCMC(nmc).L = L;
            MCMC(nmc).x = A0(:);
            MCMC(nmc).phi = B0(:);
            MCMC(nmc).theta = C0(:);
            A0 = squeeze(A_0);
            B0 = squeeze(B_0);
            C0 = squeeze(C_0);
            MCMC(nmc).a0 = A0(newvn);
            MCMC(nmc).b0 = B0(newvn);
            MCMC(nmc).c0 = C0(newvn);
            MCMC(nmc).vn = vn;
            MCMC(nmc).xv = X(vn);
            MCMC(nmc).yv = Y(vn);
            MCMC(nmc).zv = Z(vn);
            Lold = L;
            nfail2 = 0;
            siref(indoi) = si;
        else
            if ~showflag
                fid = fopen([workdir '/results/log' addname '.txt'],'at');
                fprintf(fid,'Reject current step.\n');
                fprintf(fid,'alpha = %f | crit = %f | LH = %f\n',alpha,crit,Lold);
                fclose(fid);
            else
                disp('Reject current step.')
                fprintf('alpha = %f | crit = %f | LH = %f\n',alpha,crit,Lold);
            end
            ifpmts = 1:length(X(:));
            ifpmts(vn) = [];
            noc = length(vn);
            MCMC(nmc).L = 0;
            MCMC(nmc).x = A_0(:);
            MCMC(nmc).phi = B_0(:);
            MCMC(nmc).theta = C_0(:);
            A0 = squeeze(A_0);
            B0 = squeeze(B_0);
            C0 = squeeze(C_0);
            MCMC(nmc).a0 = A0(newvn);
            MCMC(nmc).b0 = B0(newvn);
            MCMC(nmc).c0 = C0(newvn);
            MCMC(nmc).vn = vn;
            MCMC(nmc).xv = X(vn);
            MCMC(nmc).yv = Y(vn);
            MCMC(nmc).zv = Z(vn);
            nfail2 = nfail2+1;
        end
        
    end
    ll(nmc) = L;


    if showflag || printflag
        iypl = round(ny/2); 
        set(fig, 'CurrentAxes', splt1);
        if MCMC(nmc).L==0
            plot(nmc,ll(nmc),'ro');
        else
            plot(nmc,ll(nmc),'b.');
        end
        plot(nmc,L./alpha,'k.')
        plot(nmc,crit,'g.')
        plot(nmc,alpha,'r.')
        drawnow;
        apl = squeeze(A_0(:,iypl,:)).*100;
        bpl = squeeze(B_0(:,iypl,:))./pi.*180;
        cpl = squeeze(C_0(:,iypl,:))./pi.*180;
        set(fig, 'CurrentAxes', splt2);
        imagesc(x,z,apl')
        axis([min(X(:)) max(X(:)) min(Z(:)) max(Z(:))])
        ylabel('z in [km]','FontSize',9)
        if mean(apl(:)) > 0
            caxis([min(apl(:)) 100]);
        else
            caxis([-100 max(apl(:))])
        end
        cax1 = colorbar;
        ylabel(cax1,'Ratio in [pct]','FontSize',9)
        drawnow;
        set(fig, 'CurrentAxes', splt3);
        imagesc(x,z,bpl')
        axis([min(X(:)) max(X(:)) min(Z(:)) max(Z(:))])
        ylabel('z in [km]','FontSize',9)
        cax2 = colorbar;
        ylabel(cax2,'Phi in [deg]','FontSize',9)
        drawnow;
        set(fig, 'CurrentAxes', splt4);
        imagesc(x,z,cpl')
        axis([min(X(:)) max(X(:)) min(Z(:)) max(Z(:))])
        ylabel('z in [km]','FontSize',9)
        xlabel('x in [km]','FontSize',9)
        caxis([-60 60])
        cax3 = colorbar;
        ylabel(cax3,'Theta in [deg]','FontSize',9)
        drawnow;
        pause(0.1)
    end
end

if printflag
    filename2 = [filename '_' addname '_.jpg'];
    print(fig,[graphicsfolder '/' filename2],'-djpeg','-r300');
    close(fig)
    fid = fopen([workdir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Current plot printed to file.\n');
    fclose(fid);
end

Llim = max(ll).*2/3;
ind = 1:Nmc;
indoi = ind(ll>Llim);
for i = 1:length(indoi)
    aa(:,:,:,i) = reshape(MCMC(indoi(i)).x,nx,ny,nz);
    bb(:,:,:,i) = reshape(MCMC(indoi(i)).phi,nx,ny,nz);
    cc(:,:,:,i) = reshape(MCMC(indoi(i)).theta,nx,ny,nz);
    mf(i) = ll(indoi(i));
end

model.zlim = zlim;
model.edges = edgeadd;
model.dyfac = 1;
model.ix = nx;
model.iy = ny;
model.iz = nz;
model.xstep = xstep;
model.phistep = phistep;
model.thetastep = thetastep;
model.vflag = vflag;
model.thetaflag = thetaflag;
model.initflag = initflag;
model.initx = strength_in;
model.initphi = phi_in;
model.inittheta = theta_in;
model.data = data;
result.x = aa;
result.phi = bb;
result.theta = cc;
result.mf = mf;
result.MCMC = MCMC;


end