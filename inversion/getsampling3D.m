function sets = getsampling3D(data,maxfrac,vflag,dec,zlim,nx,ny,nz,edgeadd,showflag,printflag,resultsfolder,modelname)
% Copyright 2024 F.Link and M.D.Long 

% Prepare weights
X_ST=data.x.*1000;
Y_ST=data.y.*1000;
x=linspace(min(X_ST),max(X_ST),nx+1);
y=linspace(min(Y_ST),max(Y_ST),ny+1);
ind = 1:length(X_ST);
w = zeros(size(data.x));
for i = 1:nx
    for j = 1:ny
        idx = ind(X_ST>=x(i)&X_ST<=x(i+1)&Y_ST>=y(j)&Y_ST<=y(j+1));
        w(idx) = 1./length(idx);
    end
end
% compute sensitivity
[~,EK_X,~,~,EK_X2,~,~,~]=SIsensitivity(data,1:length(data.x),w,1,vflag,dec,0.1.*ones(nx,ny,nz),pi/4.*ones(nx,ny,nz),zeros(nx,ny,nz),zlim.*1000,nx,ny,nz,edgeadd,1);

xx = linspace(min(data.x).*1000-edgeadd/2,max(data.x).*1000+edgeadd/2,nx)./1000;
yy = linspace(min(data.y).*1000-edgeadd/2,max(data.y).*1000+edgeadd/2,ny)./1000;
idy = round(ny/2);
zz = linspace(min(zlim),max(zlim),nz);

if printflag
    if showflag
        fig = figure; 
    else
        fig = figure('visible','off'); 
    end
    imagesc(xx,zz,squeeze(EK_X(:,idy,:))'); 
    colorbar
    xlabel('x Profile in [km]')
    ylabel('Depth in [km]')
    title(['2D cut at y=' num2str(round(yy(idy))) ' km'])
%     print(fig,[resultsfolder '/' modelname '_sensitivity3D.jpg'],'-djpeg','-r300')
    exportgraphics(fig,[resultsfolder '/' modelname '_sensitivity3D.jpg'],'Resolution',300,'BackgroundColor','white')
    close(fig);
    if showflag
        fig = figure; 
    else
        fig = figure('visible','off');
    end
    imagesc(xx,zz,squeeze(EK_X2(:,idy,:))'); 
    colorbar
    xlabel('x Profile in [km]')
    ylabel('Depth in [km]')
    title(['2D cut at y=' num2str(round(yy(idy))) ' km'])
%     print(fig,[resultsfolder '/' modelname '_nn3D.jpg'],'-djpeg','-r300')
    exportgraphics(fig,[resultsfolder '/' modelname '_nn3D.jpg'],'Resolution',300,'BackgroundColor','white')
    close(fig)
end

% prepare grid cells
% [XX,YY,ZZ] = meshgrid(xx,yy,zz);
% XX = reshape(XX,nx,ny,nz);
% YY = reshape(YY,nx,ny,nz);
% ZZ = reshape(ZZ,nx,ny,nz);
[YY,XX,ZZ] = meshgrid(yy,xx,zz);
xn = XX(:);
yn = YY(:);
zn = ZZ(:);
ekx2 = EK_X2(:);
sets = struct('ix',cell(nx,ny,nz),'rr',cell(nx,ny,nz),'n',cell(nx,ny,nz),'ar',cell(nx,ny,nz));
rr = zeros(nx,ny,nz);
nn = zeros(nx,ny,nz);
for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            dd = sqrt((xn-xx(i)).^2+(yn-yy(j)).^2+(zn-zz(k)).^2);
            [dd,I] = sort(dd,'ascend');
            temp = cumsum(ekx2(I));
            ix = find(temp>maxfrac.*max(ekx2),1);
            ioi = I(1:ix);
            sets(i,j,k).ix = ioi;
            sets(i,j,k).rr = dd(ix);
            sets(i,j,k).n = length(ioi);
            sets(i,j,k).ar = length(ioi).*(xx(2)-xx(1)).*(yy(2)-yy(1)).*(zz(2)-zz(1));
            rr(i,j,k) = dd(ix);
            nn(i,j,k) = sets(i,j,k).n;
        end
    end
end

if printflag
    indxall = 1:length(XX(:));
    indxoi = indxall(YY(:)==yy(idy));
    indx = indxoi(round(linspace(1,length(indxoi),15)));
    if showflag
        fig = figure('Position',[50 50 1200 600]); 
    else
        fig = figure('Position',[50 50 1200 600],'visible','off'); 
    end
    imagesc(xx,zz,squeeze(nn(:,idy,:))'); 
    colorbar
    hold on
    axis equal tight
    drawnow
    ax = gca;
    AR = get(ax, 'dataaspectratio');
    if ~isequal(AR(1:2), [1 1])
      error('Units are not equal on X and Y, cannot create marker size that is one unit on both');
    end
    oldunits = get(ax, 'Units');
    set(ax, 'Units', 'points');
    pos = get(ax, 'Position');    %[X Y width height]
    set(ax, 'Units', oldunits');
    XL = xlim(ax);
    points_per_unit = pos(3) / (XL(2) - XL(1));
    marker_size = points_per_unit .^2 * pi / 4;  
    scatter(XX(indx),ZZ(indx),pi.*rr(indx).^2.*marker_size,[1 0 0])
    plot(XX(indx),ZZ(indx),'r.','MarkerSize',8)
    xlabel('x Profile in [km]')
    ylabel('Depth in [km]')
    title(['2D cut at y=' num2str(round(yy(idy))) ' km'])
%     print(fig,[resultsfolder '/' modelname '_cells3D.jpg'],'-djpeg','-r300')
    exportgraphics(fig,[resultsfolder '/' modelname '_cells3D.jpg'],'Resolution',300,'BackgroundColor','white')
    close(fig);
end


end