function [model,result] = DoInvert3Dgd(data,thetaflag,vflag,wflag,zlim,nx,ny,nz,ix,iy,iz,edgeadd,sets,xstepini,phistepini,thetastepini,N1,N2,n1,n2,nevmax,initflag,strength_in0,phi_in0,theta_in0,smoothpar,dec,showflag,printflag,graphicsfolder,filename,addname)
% Copyright 2024 F.Link and M.D.Long 

[workdir,~,~] = fileparts(graphicsfolder(1:end-1));
if ~showflag
    fid = fopen([workdir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Initialize results variables.\n');
    fclose(fid);
end

% initialize results variables
if ix~=nx || iy~=ny || iz~=nz
    aa = zeros(ix,iy,iz,N1);
    bb = aa;
    cc = aa;
else
    aa = zeros(nx,ny,nz,N1);
    bb = aa;
    cc = aa;
end
mf = zeros(length(data.x),N1);

if ~thetaflag
    n2 = 2.*N2;
end

if ~showflag
    fid = fopen([workdir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Prepare weights.\n');
    fclose(fid);
end

% Prepare weights
X_ST=data.x.*1000;
Y_ST=data.y.*1000;
x=linspace(min(X_ST),max(X_ST),nx+1);
y=linspace(min(Y_ST),max(Y_ST),ny+1);
ind = 1:length(X_ST);
if wflag
    w0 = zeros(size(data.x));
    for i = 1:nx
        for j = 1:ny
            stats(i).ind = ind(X_ST>=x(i)&X_ST<=x(i+1)&Y_ST>=y(j)&Y_ST<=y(j+1));
            w0(stats(i).ind) = 1./length(stats(i).ind);
        end
    end
else
    w0 = ones(size(data.x));
end

if ~showflag
    fid = fopen([workdir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Start loop for statistics.\n');
    fclose(fid);
end

X_ST=data.x;
Y_ST=data.y;
if ix~=nx || iy~=ny || iz~=nz
    A_0 = zeros(ix,iy,iz);
    B_0 = zeros(ix,iy,iz);
    C_0 = zeros(ix,iy,iz);
else
    A_0 = zeros(nx,ny,nz);
    B_0 = zeros(nx,ny,nz);
    C_0 = zeros(nx,ny,nz);
end
if ix~=nx || iy~=ny || iz~=nz
x=linspace(min(X_ST)-(edgeadd/1000/2),max(X_ST)+(edgeadd/1000/2),nx);
y=linspace(min(Y_ST)-(edgeadd/1000/2),max(Y_ST)+(edgeadd/1000/2),ny);
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
end


% start loop for statistics
for kk = 1:N1
    % starting model
    if ix~=nx || iy~=ny || iz~=nz
        A_0 = zeros(ix,iy,iz);
        B_0 = zeros(ix,iy,iz);
        C_0 = zeros(ix,iy,iz);
    end
    if initflag
        A_0 = strength_in0;
        B_0 = phi_in0;
        C_0 = theta_in0;
    else
        A_0(:) = 0.01;
        B_0(:) = rand.*pi;
        C_0(:) = 0;
    end
    S2 = zeros(size(data.x));
    ni = 1;
    ni2 = 1;
    
    ngflag = 1;
    if printflag
        if showflag
            fig = figure('Position',[300 50 600 580]);
        else
            fig = figure('Position',[300 50 600 580],'visible','off');
        end
    else
        fig = 0;
    end
    
    % use different subset for each statistical iteration
    noi = min(nevmax,round(length(data.x)*2/4));
    indoi = randperm(length(data.x),noi);
    aflag = 0;
    bflag = 0;
    xstep = xstepini;
    phistep = phistepini;
    thetastep = thetastepini;
    xstepo = xstep;
    phistepo = phistep;
    thetastepo = thetastep;
    xgflag = 1;
    phigflag = 1;
    if thetaflag
        thetagflag = 1;
    else
        thetagflag = 0;
    end
    
    ii = 1;
    while (xstep > 0.0075|| phistep > 0.02 || thetastep > 0.02) && ii < N2
        
        if wflag
            w = zeros(size(data.x));
            for jj = 1:length(stats)
                stats(jj).ind2 = intersect(indoi,stats(jj).ind);
                w(stats(jj).ind2) = 1./length(stats(jj).ind2);
            end
        else
            w = ones(size(data.x));
        end
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Calculate gradients.\n');
            fclose(fid);
        end
        % calculate gradients
        if ngflag
            if ix~=nx || iy~=ny || iz~=nz
                strength_in = reshape(A_0(k),nx,ny,nz);
                phi_in = reshape(B_0(k),nx,ny,nz);
                theta_in = reshape(C_0(k),nx,ny,nz);
            else
                strength_in = A_0;
                phi_in = B_0;
                theta_in = C_0;
            end
        [WF,EK_X,EK_phi,EK_theta,~,~,~,timestep]=SIsensitivity(data,indoi,w,0,vflag,dec,strength_in,phi_in,theta_in,zlim.*1000,nx,ny,nz,edgeadd,1);
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
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Done in %f seconds\n',timestep);
            fclose(fid);
        end
        ngflag = 0;
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Calculate misfit.\n');
            fclose(fid);
        end
        % and corresponding misfit
        mf2 = 1/2.*sum((L1(:)).^2+(L2(:)).^2+(L3(:)).^2)*length(WF)/4;
        misfito = 0;
        for i = 1:length(WF)
            misfito = misfito+1/2.*(WF(i).misfit.^2);
            sifin(i) = WF(i).vals;
        end
        misfito = misfito+mf2.*smoothpar.^2;
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Current misfit: %f\n',misfito);
            fclose(fid);
        else
            disp(num2str(misfito))
        end
        if printflag
            if ~showflag
                fid = fopen([workdir '/results/log' addname '.txt'],'at');
                fprintf(fid,'Prepare variables to plot.\n');
                fclose(fid);
            end
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
        if ii == 1
            silim = [0 max(abs(sioi(:)))];
            sipl0 = sipl;
            xpl0 = xpl-(xx(2)-xx(1))/10;
            ypl0 = ypl-(yy(2)-yy(1))/10;
        end
        end
        
        
        if showflag || printflag
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
        end
        end
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Step %f of %f done. (round %f)\n',ii,N2,kk);
            fclose(fid);
        else
            disp(['Step ' num2str(ii) ' of ' num2str(N2) ' done. (round ' num2str(kk) ')'])
        end
        ii = ii+1;
        if ni >= n1
            if aflag == 1
                aflag = 0;
            else
                aflag = 1;
            end
            ni = 1;
        end
        if ni2 >= n2
            if bflag == 1
                bflag = 0;
            else
                bflag = 1;
            end
            ni2 = 0;
        end

        % find model improvement
        if aflag
            if ~showflag
                fid = fopen([workdir '/results/log' addname '.txt'],'at');
                fprintf(fid,'Testing x\n');
                fclose(fid);
            else
                disp('Testing x')
            end
            EKx = transform_kernel3D(EK_X,sets);
            Rx = -(smoothpar).*max(abs(EKx(:)))./max(abs(Lx(:))).*Lx;
            Rx(isnan(Rx)) = 0; Rx(isinf(Rx)) = 0;
            EKx = EKx+Rx;
            EKx(abs(EKx)<max(abs(EKx(:)))/4) = 0;
            xfac = max(abs(EKx(:)));
            if xfac == 0
                xgflag = 0;
                xstep = xstepo;
                aflag = 0;
                ni = 1;
                if ~showflag
                    fid = fopen([workdir '/results/log' addname '.txt'],'at');
                    fprintf(fid,'No improvement in x. xstep = %f\n',xstep);
                    fclose(fid);
                else
                    disp(['No improvement in x. xstep = ' num2str(xstep)])
                end
                continue
            end
            misfit = misfito;
            nstep = 1;
            gflag = 0;
            while misfit >= misfito && nstep < 5
                ni = ni+1;
                temp = xstep.*EKx./xfac;
                if ix~=nx || iy~=ny || iz~=nz
                    temp = squeeze(temp);
                    for i = 1:length(ku)
                        temp2(i) = mean(temp(k==ku(i)));
                    end
                    temp = reshape(temp2,ix,iy,iz);
                end
                A0 = A_0+temp;
                A0(A0<0) = 0;
                A0(A0>1) = 1;
                if ix~=nx || iy~=ny || iz~=nz
                    strength_in = reshape(A0(k),nx,ny,nz);
                end
                [WF,timestep]=SIforward(data,indoi,0,vflag,strength_in,phi_in,theta_in,zlim.*1000,nx,ny,nz,edgeadd,1);
                if ~showflag
                    fid = fopen([workdir '/results/log' addname '.txt'],'at');
                    fprintf(fid,'Done in %f seconds\n',timestep);
                    fclose(fid);
                end
                sig1 = strength_in.*cos(phi_in).*cos(theta_in);
                sig2 = strength_in.*sin(phi_in).*cos(theta_in);
                sig3 = strength_in.*sin(theta_in);
                L1 = zeros(size(sig1));
                L1(2:end-1,:,2:end-1) = sig1(1:end-2,:,2:end-1)+sig1(3:end,:,2:end-1)+sig1(2:end-1,:,1:end-2)+sig1(2:end-1,:,3:end)-4.*sig1(2:end-1,:,2:end-1);
                L2 = zeros(size(sig2));
                L2(2:end-1,:,2:end-1) = sig2(1:end-2,:,2:end-1)+sig2(3:end,:,2:end-1)+sig2(2:end-1,:,1:end-2)+sig2(2:end-1,:,3:end)-4.*sig2(2:end-1,:,2:end-1);
                L3 = zeros(size(sig2));
                L3(2:end-1,:,2:end-1) = sig3(1:end-2,:,2:end-1)+sig3(3:end,:,2:end-1)+sig3(2:end-1,:,1:end-2)+sig3(2:end-1,:,3:end)-4.*sig3(2:end-1,:,2:end-1);
                mf2 = 1/2.*sum((L1(:)).^2+(L2(:)).^2+(L3(:)).^2)*length(WF)/4;
                misfit = 0;
                for i = 1:length(WF)
                    misfit = misfit+WF(i).misfit;
                end
                misfit = misfit+mf2.*smoothpar.^2;
                if ~showflag
                    fid = fopen([workdir '/results/log' addname '.txt'],'at');
                    fprintf(fid,'New misfit: %f\n',misfit);
                    fclose(fid);
                else
                    disp(num2str(misfit))
                end
                if misfit>=misfito
                    xstep = xstep/2;
                    if ~showflag
                        fid = fopen([workdir '/results/log' addname '.txt'],'at');
                        fprintf(fid,'xstep = %f\n',xstep);
                        fclose(fid);
                    else
                        disp(['xstep = ' num2str(xstep)])
                    end
                else
                    ngflag = 1;
                    xgflag = 1;
                    misfito = misfit;
                    xstepo = xstep;
                    A_0 = A0;
                    gflag = 1;
                    break
                end
                nstep = nstep+1;
            end
            if ~gflag
                xgflag = 0;
                xstep = xstepo;
                aflag = 0;
                ni = 1;
                if ~showflag
                    fid = fopen([workdir '/results/log' addname '.txt'],'at');
                    fprintf(fid,'No improvement in x. xstep = %f\n',xstep);
                    fclose(fid);
                else
                    disp(['No improvement in x. xstep = ' num2str(xstep)])
                end
            end
        else
            if ~showflag
                fid = fopen([workdir '/results/log' addname '.txt'],'at');
                fprintf(fid,'Testing phi\n');
                fclose(fid);
            else
                disp('Testing phi')
            end
            EKphi = transform_kernel3D(EK_phi,sets);
            Rp = -(smoothpar).*max(abs(EKphi(:)))./max(abs(Lp(:))).*Lp;
            Rp(isnan(Rp)) = 0; Rp(isinf(Rp)) = 0;
            EKphi = EKphi+Rp;
            phifac = max(abs(EKphi(:)));
            if phifac == 0
                phigflag = 0;
                phistep = phistepo;
                aflag = 1;
                ni = 1;
                if ~showflag
                    fid = fopen([workdir '/results/log' addname '.txt'],'at');
                    fprintf(fid,'No improvement in phi. phistep = %f\n',phistep);
                    fclose(fid);
                else
                    disp(['No improvement in phi. phistep = ' num2str(phistep)])
                end
                continue
            end
            misfit = misfito;
            nstep = 1;
            gflag = 0;
            while misfit >= misfito && nstep < 5
                ni = ni+1;
                temp = phistep.*EKphi./phifac;
                if ix~=nx || iy~=ny || iz~=nz
                    temp = squeeze(temp);
                    for i = 1:length(ku)
                        temp2(i) = mean(temp(k==ku(i)));
                    end
                    temp = reshape(temp2,ix,iy,iz);
                end
                B0 = mod(B_0+temp,pi);
                if ix~=nx || iy~=ny || iz~=nz
                    phi_in = reshape(B0(k),nx,ny,nz);
                end
                [WF,timestep]=SIforward(data,indoi,0,vflag,strength_in,phi_in,theta_in,zlim.*1000,nx,ny,nz,edgeadd,1);
                if ~showflag
                    fid = fopen([workdir '/results/log' addname '.txt'],'at');
                    fprintf(fid,'Done in %f seconds\n',timestep);
                    fclose(fid);
                end
                sig1 = strength_in.*cos(phi_in).*cos(theta_in);
                sig2 = strength_in.*sin(phi_in).*cos(theta_in);
                sig3 = strength_in.*sin(theta_in);
                L1 = zeros(size(sig1));
                L1(2:end-1,:,2:end-1) = sig1(1:end-2,:,2:end-1)+sig1(3:end,:,2:end-1)+sig1(2:end-1,:,1:end-2)+sig1(2:end-1,:,3:end)-4.*sig1(2:end-1,:,2:end-1);
                L2 = zeros(size(sig2));
                L2(2:end-1,:,2:end-1) = sig2(1:end-2,:,2:end-1)+sig2(3:end,:,2:end-1)+sig2(2:end-1,:,1:end-2)+sig2(2:end-1,:,3:end)-4.*sig2(2:end-1,:,2:end-1);
                L3 = zeros(size(sig2));
                L3(2:end-1,:,2:end-1) = sig3(1:end-2,:,2:end-1)+sig3(3:end,:,2:end-1)+sig3(2:end-1,:,1:end-2)+sig3(2:end-1,:,3:end)-4.*sig3(2:end-1,:,2:end-1);
                mf2 = 1/2.*sum((L1(:)).^2+(L2(:)).^2+(L3(:)).^2)*length(WF)/4;
                misfit = 0;
                for i = 1:length(WF)
                    misfit = misfit+WF(i).misfit;
                end
                misfit = misfit+mf2.*smoothpar.^2;
                if ~showflag
                    fid = fopen([workdir '/results/log' addname '.txt'],'at');
                    fprintf(fid,'New misfit: %f\n',misfit);
                    fclose(fid);
                else
                    disp(num2str(misfit))
                end
                if misfit>=misfito
                    phistep = phistep/2;
                    if ~showflag
                        fid = fopen([workdir '/results/log' addname '.txt'],'at');
                        fprintf(fid,'phistep = %f\n',phistep);
                        fclose(fid);
                    else
                        disp(['phistep = ' num2str(phistep)])
                    end
                else
                    ngflag = 1;
                    phigflag = 1;
                    misfito = misfit;
                    phistepo = phistep;
                    B_0 = B0;
                    gflag = 1;
                    break
                end
                nstep = nstep+1;
            end
            if ~gflag
                phigflag = 0;
                phistep = phistepo;
                aflag = 1;
                ni = 1;
                if ~showflag
                    fid = fopen([workdir '/results/log' addname '.txt'],'at');
                    fprintf(fid,'No improvement in phi. phistep = %f\n',phistep);
                    fclose(fid);
                else
                    disp(['No improvement in phi. phistep = ' num2str(phistep)])
                end
            end
        end
        if bflag
            if ~showflag
                fid = fopen([workdir '/results/log' addname '.txt'],'at');
                fprintf(fid,'Testing theta\n');
                fclose(fid);
            else
                disp('Testing theta')
            end
            EKtheta = transform_kernel3D(EK_theta,sets);
            Rt = -(smoothpar).*max(abs(EKtheta(:)))./max(abs(Lt(:))).*Lt;
            Rt(isnan(Rt)) = 0; Rt(isinf(Rt)) = 0;
            EKtheta = EKtheta+Rt;
            EKtheta(abs(EKtheta)<max(abs(EKtheta(:)))/4*2) = 0;
            thetafac = max(abs(EKtheta(:)));
            if thetafac == 0
                thetagflag = 0;
                thetastep = thetastepo;
                bflag = 0;
                ni2 = 1;
                if ~showflag
                    fid = fopen([workdir '/results/log' addname '.txt'],'at');
                    fprintf(fid,'No improvement in theta. thetastep = %f\n',thetastep);
                    fclose(fid);
                else
                    disp(['No improvement in theta. thetastep = ' num2str(thetastep)])
                end
                continue
            end
            misfit = misfito;
            nstep = 1;
            gflag = 0;
            while misfit >= misfito && nstep < 5
                temp = thetastep.*EKtheta./thetafac;
                if ix~=nx || iy~=ny || iz~=nz
                    temp = squeeze(temp);
                    for i = 1:length(ku)
                        temp2(i) = mean(temp(k==ku(i)));
                    end
                    temp = reshape(temp2,ix,iy,iz);
                end
                C0 = C_0-temp;
                C0(C0<-pi/3) = -pi/3;
                C0(C0>pi/3) = pi/3;
                if ix~=nx || iz~=nz
                    theta_in = reshape(C0(k),nx,ny,nz);
                end
                [WF,timestep]=SIforward(data,indoi,0,vflag,strength_in,phi_in,theta_in,zlim.*1000,nx,ny,nz,edgeadd,1);
                if ~showflag
                    fid = fopen([workdir '/results/log' addname '.txt'],'at');
                    fprintf(fid,'Done in %f seconds\n',timestep);
                    fclose(fid);
                end
                sig1 = strength_in.*cos(phi_in).*cos(theta_in);
                sig2 = strength_in.*sin(phi_in).*cos(theta_in);
                sig3 = strength_in.*sin(theta_in);
                L1 = zeros(size(sig1));
                L1(2:end-1,:,2:end-1) = sig1(1:end-2,:,2:end-1)+sig1(3:end,:,2:end-1)+sig1(2:end-1,:,1:end-2)+sig1(2:end-1,:,3:end)-4.*sig1(2:end-1,:,2:end-1);
                L2 = zeros(size(sig2));
                L2(2:end-1,:,2:end-1) = sig2(1:end-2,:,2:end-1)+sig2(3:end,:,2:end-1)+sig2(2:end-1,:,1:end-2)+sig2(2:end-1,:,3:end)-4.*sig2(2:end-1,:,2:end-1);
                L3 = zeros(size(sig2));
                L3(2:end-1,:,2:end-1) = sig3(1:end-2,:,2:end-1)+sig3(3:end,:,2:end-1)+sig3(2:end-1,:,1:end-2)+sig3(2:end-1,:,3:end)-4.*sig3(2:end-1,:,2:end-1);
                mf2 = 1/2.*sum((L1(:)).^2+(L2(:)).^2+(L3(:)).^2)*length(WF)/4;
                misfit = 0;
                for i = 1:length(WF)
                    misfit = misfit+WF(i).misfit;
                end
                misfit = misfit+mf2.*smoothpar.^2;
                if ~showflag
                    fid = fopen([workdir '/results/log' addname '.txt'],'at');
                    fprintf(fid,'New misfit: %f\n',misfit);
                    fclose(fid);
                else
                    disp(num2str(misfit))
                end
                if misfit>=misfito
                    thetastep = thetastep/2;
                    if ~showflag
                        fid = fopen([workdir '/results/log' addname '.txt'],'at');
                        fprintf(fid,'thetastep = %f\n',thetastep);
                        fclose(fid);
                    else
                        disp(['thetastep = ' num2str(thetastep)])
                    end
                else
                    ngflag = 1;
                    thetagflag = 1;
                    misfito = misfit;
                    thetastepo = thetastep;
                    C_0 = C0;
                    gflag = 1;
                    break
                end
                nstep = nstep+1;
            end
            if ~gflag
                thetagflag = 0;
                thetastep = thetastepo;
                bflag = 0;
                ni2 = 1;
                if ~showflag
                    fid = fopen([workdir '/results/log' addname '.txt'],'at');
                    fprintf(fid,'No improvement in theta. thetastep = %f\n',thetastep);
                    fclose(fid);
                else
                    disp(['No improvement in theta. thetastep = ' num2str(thetastep)])
                end
            end
        end
        ni2 = ni2+1;
        if ~xgflag && ~phigflag && ~thetagflag
            break
        end
    end
    eflag(:,kk) = zeros(size(data.x));
    eflag(indoi,kk) = 1;
    S2(indoi) = sifin;
    aa(:,:,:,kk) = A_0;
    bb(:,:,:,kk) = mod(B_0,pi);
    cc(:,:,:,kk) = C_0;
    mf(:,kk) = S2;
    if ~showflag
        fid = fopen([workdir '/results/log' addname '.txt'],'at');
        fprintf(fid,'Current step stored in variables.\n');
        fclose(fid);
    end
    if printflag
        filename2 = [filename '_' addname num2str(kk) '_.jpg'];
        print(fig,[graphicsfolder '/' filename2],'-djpeg','-r300');
        close(fig)
        fid = fopen([workdir '/results/log' addname '.txt'],'at');
        fprintf(fid,'Current plot printed to file.\n');
        fclose(fid);
    end
end

model.sets = sets;
model.zlim = zlim;
model.edges = edgeadd;
model.dyfac = 1;
model.nx = nx;
model.ny = ny;
model.nz = nz;
model.ix = ix;
model.iy = iy;
model.iz = iz;
model.xstep = xstep;
model.phistep = phistep;
model.thetastep = thetastep;
model.N1 = N1;
model.N2 = N2;
model.n1 = n1;
model.n2 = n2;
model.vflag = vflag;
model.thetaflag = thetaflag;
model.initflag = initflag;
model.initx = strength_in0;
model.initphi = phi_in0;
model.inittheta = theta_in0;
model.data = data;
model.w = w0;
result.eflag = eflag;
result.x = aa;
result.phi = bb;
result.theta = cc;
result.mf = mf;

end