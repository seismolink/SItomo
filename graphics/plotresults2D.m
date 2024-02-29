function plotresults2D(sel_data)
% Copyright 2024 F.Link and M.D.Long 

if ~isstruct(sel_data)
    input_flag = 0;
    clear sel_data
else
    if ~sel_data.remoteflag
        input_flag = 0;
    else
        input_flag = 1;
    end
end

% Set FontSize for Axes
AFS = 9;
% Set FontSize for Titles
TFS = 11;
% Set FontSize for Legend
LFS = 8;

if ~input_flag
    sel_data.SI_dir = input('Define working directory for SI tomography: ','s');
end

if ~input_flag
    minval = input('Define minimum value for X in the final model: ');
else
    minval = sel_data.minval;
end

col1 = copper(64);
col2 = [linspace(col1(end,1),1,8)' linspace(col1(end,2),1,8)' linspace(col1(end,3),1,8)'];
colp = [col1; col2;];

if ~exist([sel_data.SI_dir '/input/input2D.mat'],'file')
    disp('No input structure. Import splitting intensity data first for 2D structure.')
    return
else
    load([sel_data.SI_dir '/input/input2D.mat']);
end

if exist([sel_data.SI_dir '/results/finmod.mat'],'file')

    % load weighted mean result
    load([sel_data.SI_dir '/results/finmod.mat'])
    xoi = finmod.x;
    ex = finmod.ex;
    phif = finmod.phi;
    ephif = finmod.ephi;
    thetaf = finmod.theta;
    ethetaf = finmod.etheta;
    esi = finmod.esi;
    si = finmod.si;
    data = finmod.data;
    model = finmod.model;
    xpl = finmod.xx;
    zpl = finmod.zz
    indx = 1:length(si);
    clear finmod
    zlim = model.zlim;
    
else
    
    A = dir([sel_data.SI_dir '/results/result*']);
    stflag = 1;
    for i = 1:length(A)
        if contains(A(i).name,'result3D')
            continue
        end
        load([A(i).folder '/' A(i).name]);
        fn = fieldnames(result);
        if sum(strcmp(fn,'MCMC'))>0
            mcmcflag = 1;
            clear lmf
            for j = 1:length(result.MCMC)-1
                lmf(j) = result.MCMC(j).L;
            end
            [~,I] = sort(lmf,'descend');
            len = sum(lmf~=0);
            disp(num2str(len))
            fr = 5;
            if stflag == 1
                M = size(result.x);
                mf = lmf(I(1:round(len/fr)));
                all.mf = mf;
                n = 1;
                m = 1;
                for j = 1:length(I(1:round(len/fr)))
                    N = min(length(result.MCMC(I(j)).vn),length(result.MCMC(I(j)).a0));
                    if length(result.MCMC(I(j)).vn) ~= length(result.MCMC(I(j)).a0)
                        todel(m) = j;
                        m = m+1;
                        continue
                    end
                    all.vn(:,n) = nan(1000,1);
                    all.xn(:,n) = nan(1000,1);
                    all.zn(:,n) = nan(1000,1);
                    all.a0(:,n) = nan(1000,1);
                    all.b0(:,n) = nan(1000,1);
                    all.c0(:,n) = nan(1000,1);
                    all.vn(1:N,n) = result.MCMC(I(j)).vn(1:N);
                    all.xn(1:N,n) = result.MCMC(I(j)).xv(1:N);
                    all.zn(1:N,n) = result.MCMC(I(j)).zv(1:N);
                    all.a0(1:N,n) = result.MCMC(I(j)).a0(1:N);
                    all.b0(1:N,n) = result.MCMC(I(j)).b0(1:N);
                    all.c0(1:N,n) = result.MCMC(I(j)).c0(1:N);
                    n = n+1;
                end
                all.mf(todel) = [];
                stflag = 0;
            else
                clear todel
                M2 = size(all.mf);
                all.mf(M2(end)+1:M2(end)+round(len/fr)) = lmf(I(1:round(len/fr)));
                n = 1;
                m = 1;
                for j = 1:length(I(1:round(len/fr)))
                    N = min(length(result.MCMC(I(j)).vn),length(result.MCMC(I(j)).a0));
                    if length(result.MCMC(I(j)).vn) ~= length(result.MCMC(I(j)).a0)
                        todel(m) = M2(end)+j;
                        m = m+1;
                        continue
                    end
                    all.vn(:,n) = nan(1000,1);
                    all.xn(:,n) = nan(1000,1);
                    all.zn(:,n) = nan(1000,1);
                    all.a0(:,n) = nan(1000,1);
                    all.b0(:,n) = nan(1000,1);
                    all.c0(:,n) = nan(1000,1);
                    all.vn(1:N,M2(end)+n) = result.MCMC(I(j)).vn(1:N);
                    all.xn(1:N,M2(end)+n) = result.MCMC(I(j)).xv(1:N);
                    all.zn(1:N,M2(end)+n) = result.MCMC(I(j)).zv(1:N);
                    all.a0(1:N,M2(end)+n) = result.MCMC(I(j)).a0(1:N);
                    all.b0(1:N,M2(end)+n) = result.MCMC(I(j)).b0(1:N);
                    all.c0(1:N,M2(end)+n) = result.MCMC(I(j)).c0(1:N);
                    n = n+1;
                end
                all.mf(todel) = [];
            end
        else
            mcmcflag = 0;
            if stflag == 1
                all.x = result.x;
                all.phi = result.phi;
                all.theta = result.theta;
                [MM,NN] = size(result.mf);
                if MM == length(model.data.si)
                    for k = 1:NN
                        ioi = abs(result.mf(:,k))>eps;
                        mf(k) = sum(0.5.*((model.data.si(ioi)'-result.mf(ioi,k)).^2));
                    end
                    mf = mf';
                else
                    mf = result.mf;
                end
                all.mf = mf;
                stflag = 0; 
            else
                M = size(result.x);
                M2 = size(all.x);
                [MM,NN] = size(result.mf);
                if length(M) == 2
                    M(3) = 1;
                end
                if length(M2) == 2
                    M2(3) = 1;
                end
                all.x(:,:,M2(end)+1:M2(end)+M(end)) = result.x;
                all.phi(:,:,M2(end)+1:M2(end)+M(end)) = result.phi;
                all.theta(:,:,M2(end)+1:M2(end)+M(end)) = result.theta;
                if MM == length(model.data.si)
                    for k = 1:NN
                        ioi = result.mf(:,k)~=0;
                        mf(k) = sum(0.5.*((model.data.si(ioi)'-result.mf(ioi,k)).^2));
                    end
                    mf = mf';
                else
                    mf = result.mf;
                end
                all.mf(M2(end)+1:M2(end)+M(end),:) = mf;

            end
        end
    end
    result = all;
    if length(result.mf) > 1
    [~,I] = sort(result.mf,'ascend');
    gflag = abs(result.mf)<=abs(result.mf(I(round(length(result.mf)./5.*4))));
    weights = 1-(abs(result.mf)./(max(abs(result.mf))).*1.5);
    weights(~gflag) = 0;
    
    if ~mcmcflag
        result.x(:,:,weights==0) = [];
        result.phi(:,:,weights==0) = [];
        result.theta(:,:,weights==0) = [];
        result.mf(weights==0) = [];
    else
        result.vn(:,weights==0) = [];
        result.xn(:,weights==0) = [];
        result.zn(:,weights==0) = [];
        result.a0(:,weights==0) = [];
        result.b0(:,weights==0) = [];
        result.c0(:,weights==0) = [];
        xx = result.xn(:); xx(isnan(xx)) = [];
        zz = result.zn(:); zz(isnan(zz)) = [];
        a0 = result.a0(:); a0(isnan(a0)) = [];
        b0 = result.b0(:); b0(isnan(b0)) = [];
        c0 = result.c0(:); c0(isnan(c0)) = [];
        todel = ((abs(xx)+abs(zz)+abs(a0)+abs(b0)+abs(c0)==0)+a0>0.5)>0;
        xx(todel) = [];
        zz(todel) = [];
        a0(todel) = [];
        b0(todel) = [];
        c0(todel) = [];
        x = linspace(min(xx),max(xx),model.ix);
        z = linspace(min(zz),max(zz),model.iz);
        ddx = x(2)-x(1);
        [X,Z] = meshgrid(x,z);
        xn = X(:);
        zn = Z(:);
        for i = 1:length(xn)
            dd(i,:) = sqrt((xx-xn(i)).^2+(zz-zn(i)).^2);
        end
        w = 1./dd.^3;
        w(dd>(3*ddx)) = 0;
        wu = unique(w);
        w(isinf(w)) = wu(end-1);
        for i = 1:1000
            ioi = randi(length(b0),round((M(3).*M(2).*M(1))./2),1);
            A_0(:,:,i) = reshape(a0(ioi)'*w(:,ioi)'./sum(w(:,ioi),2)',model.iz,model.ix)';
            b0a = reshape(sin(2*b0(ioi))'*w(:,ioi)'./sum(w(:,ioi),2)',model.iz,model.ix)';
            b0b = reshape(cos(2*b0(ioi))'*w(:,ioi)'./sum(w(:,ioi),2)',model.iz,model.ix)';
            B_0(:,:,i) = atan2(b0a,b0b)/2;
            C_0(:,:,i) = reshape(c0(ioi)'*w(:,ioi)'./sum(w(:,ioi),2)',model.iz,model.ix)';
        end
        A_0(isnan(A_0)) = 0;
        B_0(isnan(B_0)) = 0;
        C_0(isnan(C_0)) = 0;
        result.x = A_0;%reshape(A_0,M(1),M(2),M(3),1000);
        result.phi = B_0;%reshape(B_0,M(1),M(2),M(3),1000);
        result.theta = C_0;%reshape(C_0,M(1),M(2),M(3),1000);
    end
    end
    in = result.x<0;
    result.x(in) = abs(result.x(in));
    result.phi(in) = mod(result.phi(in)+pi/2,pi);

    M = size(result.x);
    xn = result.x;
    phin = result.phi;
    thetan = result.theta;
    xn2 = xn;
    phin2 = phin;
    thetan2 = thetan;

    if length(M)>2
        % correct for outliers
        xedges = linspace(0,1,100);
        phiedges = linspace(0,pi,180);
        thetaedges = linspace(-pi/2,pi/2,180);
        nn = 2;
        for k = 1:M(3)
            for i = 1:M(1)
                for j = 1:M(2)
                    if xn(i,j,k) == 0
                        xn2(i,j,k) = 0;
                        phin2(i,j,k) = 0;
                        thetan2(i,j,k) = 0;
                    else
                    i1 = max(1,i-nn);
                    i2 = min(M(1),i+nn);
                    i3 = max(1,j-nn);
                    i4 = min(M(2),j+nn);
                    temp = xn(i1:i2,i3:i4,k);
                    temp = temp(:);
                    indx = temp==0;
                    temp(indx) = [];
                    [N,X] = histcounts(temp,xedges);
                    [~,I] = max(N);
                    xn2(i,j,k) = mean(X(I:I+1));
                    xn2(i,j,k) = mean(temp);
                    temp = phin(i1:i2,i3:i4,k);
                    temp = temp(:);
                    temp(indx) = [];
                    ephi1 = std(temp);
                    temp2 = temp;
                    temp2(temp>pi/2) = temp2(temp>pi/2)-pi;
                    ephi2 = std(temp2);
                    if ephi2<ephi1
                        temp = temp2;
                    end
                    [N,Phi] = histcounts(temp,phiedges);
                    [~,I] = max(N);
                    phin2(i,j,k) = mod(mean(Phi(I:I+1)),pi);
                    phin2(i,j,k) = mod(mean(temp),pi);
                    temp = thetan(i1:i2,i3:i4,k);
                    temp = temp(:);
                    temp(indx) = [];
                    [N,Theta] = histcounts(temp,thetaedges);
                    [~,I] = max(N);
                    thetan2(i,j,k) = mean(Theta(I:I+1));
                    thetan2(i,j,k) = mean(temp);
                    end
                end
            end
        end
        xx = mean(xn2,3);
        ex = std(xn2,[],3);
        xoi = xx-ex;
        xoi(xoi<minval) = 0;
        xoi(xoi~=0) = xx(xoi~=0);
        pphi = sum(phin2.*xn2,3)./sum(xn2,3);
        ephi = sqrt(sum(xn2.*(phin2-pphi).^2,3)./((sum(xn2>0.01,3)-1)./sum(xn2>0.01,3).*sum(xn2,3)));
        ttheta = sum(thetan2.*xn2,3)./sum(xn2,3);
        etheta = sqrt(sum(xn2.*(thetan2-ttheta).^2,3)./((sum(xn2>0.01,3)-1)./sum(xn2>0.01,3).*sum(xn2,3)));

        phis = mod(result.phi,pi);
        phis(phis>pi/2) = phis(phis>pi/2)-pi;
        phin = sum(phis.*xn2,3)./sum(xn2,3);
        ephin = sqrt(sum(xn2.*(phis-phin).^2,3)./((sum(xn2>0.01,3)-1)./sum(xn2>0.01,3).*sum(xn2,3)));
        phif = pphi;
        phif(ephi>ephin) = phin(ephi>ephin);
        phif = mod(phif,pi);
        ephif = ephi;
        ephif(ephi>ephin) = ephin(ephi>ephin);
        phif(xoi==0) = 0;
        phifo = phif;
        phif = phifo-(90-inp(1).caz)./180.*pi;
        phif(xoi==0) = 0;
        thetaf = ttheta;
        thetaf(xoi==0) = 0;
        ethetaf = etheta;
    else
        xoi = result.x;
        xoi(xoi<minval) = 0;
        ex = zeros(size(xoi));
        phifo = result.phi;
        phifo(xoi==0) = 0;
        ephif = ex;
        phif = mod(phifo-(90-inp(1).caz)./180.*pi,pi);
        phif(xoi==0) = 0;
        thetaf = result.theta;
        thetaf(xoi==0) = 0;
        ethetaf = ex;
    end
    
    if ~mcmcflag
    nx = model.nx;
    ny = model.ny;
    nz = model.nz;
    ix = model.ix;
    iz = model.iz;
    nny = model.dyfac;
    zlim = model.zlim;
    edgeadd = model.edges;
    data = model.data;
    todel = isnan(data.si);
    data.x(todel) = [];
    data.y(todel) = [];
    data.baz(todel) = [];
    data.az(todel) = [];
    data.per(todel) = [];
    data.p(todel) = [];
    data.si(todel) = [];
    if nx ~= ix || nz~=iz
        X_ST = data.x;
        x=linspace(min(X_ST)-(edgeadd/1000/2),max(X_ST)+(edgeadd/1000/2),nx);
        z=linspace(zlim(1),zlim(2),nz);
        xn=linspace(min(x),max(x),ix+1);
        zn=linspace(min(z),max(z),iz+1);
        dxn = xn(2)-xn(1);
        dzn = zn(2)-zn(1);
        xn(end) = [];
        zn(end) = [];
        xn = xn+dxn/2;
        zn = zn+dzn/2;
        [X,Z] = meshgrid(x,z);
        [Xn,Zn] = meshgrid(xn,zn);
        Xn = Xn'; Zn = Zn';X = X'; Z = Z';
        k = dsearchn([Xn(:) Zn(:)],[X(:) Z(:)]);
        A0 = reshape(xoi(k),nx,ny,nz);
        B0 = reshape(phifo(k),nx,ny,nz);
        C0 = reshape(thetaf(k),nx,ny,nz);
        xpl = xn;
        zpl = zn;
    else
        A0 = reshape(xoi,nx,ny,nz);
        B0 = reshape(phifo,nx,ny,nz);
        C0 = reshape(thetaf,nx,ny,nz);
        X_ST = data.x;
        xpl=linspace(min(X_ST)-(edgeadd/1000/2),max(X_ST)+(edgeadd/1000/2),nx);
        zpl=linspace(zlim(1),zlim(2),nz);
    end
    else
        A0 = reshape(xoi,model.ix,1,model.iz);
        B0 = reshape(phifo,model.ix,1,model.iz);
        C0 = reshape(thetaf,model.ix,1,model.iz);
        data = model.data;
        X_ST = data.x;
        edgeadd = model.edges;
        zlim = model.zlim;
        nx = model.ix;
        ny = model.iy;
        nz = model.iz;
        nny = model.dyfac;
        xpl=linspace(min(X_ST)-(edgeadd/1000/2),max(X_ST)+(edgeadd/1000/2),model.ix);
        zpl=linspace(zlim(1),zlim(2),model.iz);
    end
    idx = 1:length(data.si);
    [WF]=SIforward(data,idx,0,0,A0,B0,C0,zlim.*1000,nx,ny,nz,edgeadd,nny);
    for i = 1:length(WF); si(i) = WF(i).vals; end

    esi = data.si(idx)-si;

    % save weighted mean result
    finmod.x = xoi;
    finmod.ex = ex;
    finmod.phi = phif;
    finmod.ephi = ephif;
    finmod.theta = thetaf;
    finmod.etheta = ethetaf;
    finmod.phiorig = phifo;
    finmod.esi = esi;
    finmod.si = si;
    finmod.data = data;
    finmod.model = model;
    finmod.xx = xpl;
    finmod.zz = zpl;
    save([sel_data.SI_dir '/results/finmod.mat'],'finmod');
    clear finmod
end

xred = 0;
x = xpl;
z = zpl;

fig1 = figure; 
imagesc(x,z,mod(squeeze(phif)'/pi*180,180)); 
colormap(hsv);
cb = colorbar; 
caxis([0 180])
ylabel(cb,'Fast axis in [deg]','FontSize',AFS);
axis equal; 
axis([min(data.x)+xred max(data.x)-xred zlim(1) zlim(2)]); 
ylabel('z in [km]','FontSize',AFS)
xlabel('x in [km]','FontSize',AFS)
print(fig1,[sel_data.SI_dir '/graphics/result_Phi.jpg'],'-djpeg','-r600')
close(fig1)

fig2 = figure; 
imagesc(x,z,(squeeze(xoi)').*100); 
colormap(flip(colp));
cb = colorbar;
caxis([0 ceil(max(xoi(:)).*100)])
ylabel(cb,'Strength of Anisotropy in [pct]','FontSize',AFS);
axis equal; 
axis([min(data.x)+xred max(data.x)-xred zlim(1) zlim(2)]); 
ylabel('z in [km]','FontSize',AFS)
xlabel('x in [km]','FontSize',AFS)
print(fig2,[sel_data.SI_dir '/graphics/result_X.jpg'],'-djpeg','-r600')
close(fig2)

fig2 = figure; 
imagesc(x,z,squeeze(thetaf)'/pi*180); 
colormap(jet);
cb = colorbar;
ylabel(cb,'Plunging Axis in [deg]','FontSize',AFS);
caxis([-60 60])
axis equal; 
axis([min(data.x)+xred max(data.x)-xred zlim(1) zlim(2)]); 
ylabel('z in [km]','FontSize',AFS)
xlabel('x in [km]','FontSize',AFS)
print(fig2,[sel_data.SI_dir '/graphics/result_Theta.jpg'],'-djpeg','-r600')
close(fig2)


fig3 = figure; 
[~,edges] = histcounts(data.si(idx),50);
histogram(esi,'BinEdges',edges); 
hold on; 
histogram(data.si(idx),'BinEdges',edges)
legend('Residuals','Splitting Intensities','FontSize',LFS)
xlabel('Splitting Intensity Error','FontSize',AFS)
print(fig3,[sel_data.SI_dir '/graphics/SI_Variation.jpg'],'-djpeg','-r600')
close(fig3)


fig4 = figure; 
plot(data.si(idx),'.'); 
hold on; 
plot(si,'ro')
axis([0 length(data.si(idx))+1 -(max(max(abs(data.si(idx))),max(abs(si)))+0.1) max(max(abs(data.si(idx))),max(abs(si)))+0.1]) 
xlabel('Event count (sorted for receiver location','FontSize',AFS)
ylabel('Splitting Intensity','FontSize',AFS)
print(fig4,[sel_data.SI_dir '/graphics/result_SI.jpg'],'-djpeg','-r600')
close(fig4)

fig5 = figure; 
plot(esi,'.')
axis([0 length(data.si(idx))+1 -(max(max(abs(data.si(idx))),max(abs(si)))+0.1) max(max(abs(data.si(idx))),max(abs(si)))+0.1])
xlabel('Event count (sorted for receiver location','FontSize',AFS)
ylabel('SI Deviation','FontSize',AFS)
print(fig5,[sel_data.SI_dir '/graphics/result_SI_corrected.jpg'],'-djpeg','-r600')
close(fig5)

fig10 = figure; 
imagesc(x,z,mod(squeeze(ephif)'/pi*180,180)); 
colormap(jet);
cb = colorbar; 
caxis([0 90])
ylabel(cb,'Error in Fast Axis [deg]','FontSize',AFS);
axis equal; 
axis([min(data.x)+xred max(data.x)-xred zlim(1) zlim(2)]); 
ylabel('z in [km]','FontSize',AFS)
xlabel('x in [km]','FontSize',AFS)
print(fig10,[sel_data.SI_dir '/graphics/error_Phi.jpg'],'-djpeg','-r600')
close(fig10)

fig11 = figure; 
imagesc(x,z,(squeeze(ex)').*100); 
colormap(flip(colp));
cb = colorbar; 
caxis([0 max(xoi(:))])
ylabel(cb,'Error in Anisotropy [pct]','FontSize',AFS);
axis equal; 
axis([min(data.x)+xred max(data.x)-xred zlim(1) zlim(2)]); 
ylabel('z in [km]','FontSize',AFS)
xlabel('x in [km]','FontSize',AFS)
print(fig11,[sel_data.SI_dir '/graphics/error_X.jpg'],'-djpeg','-r600')
close(fig11)

fig11 = figure; 
imagesc(x,z,squeeze(ethetaf)'/pi*180); 
colormap(jet);
cb = colorbar; 
caxis([0 60])
ylabel(cb,'Error in Plunge [deg]','FontSize',AFS);
axis equal; 
axis([min(data.x)+xred max(data.x)-xred zlim(1) zlim(2)]); 
ylabel('z in [km]','FontSize',AFS)
xlabel('x in [km]','FontSize',AFS)
print(fig11,[sel_data.SI_dir '/graphics/error_Theta.jpg'],'-djpeg','-r600')
close(fig11)

end