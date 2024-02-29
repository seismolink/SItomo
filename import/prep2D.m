function prep2D(sel_data)
% prepare data for 2D inversion

% Copyright 2024 F.Link and M.D.Long 

if ~isstruct(sel_data) 
    input_flag = 0;
    clear sel_data
else
    input_flag = 1;
end

psflag = 0;

if ~input_flag || ~isfield(sel_data,'SI_dir')
    sel_data.SI_dir = input('Define working directory for SI tomography: ','s');
    dmax = input('Maximum distance of station to profile in degree (e.g., 0.6): ');
end


load([sel_data.SI_dir '/input/SIdata.mat']);
if ~exist('data','var')
    data = RES;
end
n = 1;
m = 1;
for i = 1:length(data)
    if isempty(data(i).lat)
        todel(m) = i;
        m = m+1;
        continue
    end
    for j = 1:length(data(i).ev)
        lat(n) = data(i).lat(1);
        lon(n) = data(i).lon(1);
        n = n+1;
    end
end
if m>1
    data(todel) = [];
end

idx = 1:length(lat);

minlat = min(lat);
minlon = min(lon);
maxlat = max(lat);
maxlon = max(lon);
rlat = (minlat-dmax):(dmax):(maxlat+dmax);
rlon = (minlon-dmax):(dmax):(maxlon+dmax);
plat = ((minlat-dmax/2):(dmax):(maxlat+dmax/2));
plon = ((minlon-dmax/2):(dmax):(maxlon+dmax/2));
rho = zeros(length(rlat)-1,length(rlon)-1);
for i = 1:length(rlat)-1
    for j = 1:length(rlon)-1
        nn = sum(lat<=rlat(i+1)&lat>=rlat(i)&lon<=rlon(j+1)&lon>=rlon(j));
        rho(i,j) = nn;
    end
end
rho(rho==0) = nan;

figure;
% s = pcolor(plon,plat,rho); 
% s.EdgeColor = 'none';
s = imagesc(plon,plat,rho);
set(gca,'Ydir','normal')
colormap(jet)
cb = colorbar;
ylabel(cb,'Number of Events')
hold on
plot(lon,lat,'w.','MarkerSize',10)
plot(lon,lat,'k.')
axis equal;
axis([min(rlon) max(rlon) min(rlat) max(rlat)]); 
hold on; 
grid on
xlabel('Longitude');
ylabel('Latitude');
fininp = 0;
lat1 = nan;
lon1 = nan;
lat2 = nan;
lon2 = nan;
while ~strcmp(fininp,'y')
    fininp = input('Define start point (s), end point (e) or finish selection (y): ','s');
    if strcmp(fininp,'s')
        if exist('pl1','var')
            delete(pl1);
            delete(pl2);
        end
        [x,y] = ginput(1);
        dd = sqrt((lon-x).^2+(lat-y).^2);
        [~,Imin] = min(dd);
        ddmin = sqrt((lat(Imin)-lat(idx)).^2+(lon(Imin)-lon(idx)).^2);
        lat1 = mean(lat(idx(ddmin<0.5)));
        lon1 = mean(lon(idx(ddmin<0.5)));
        pl1 = plot(lon(Imin),lat(Imin),'ro');
        pl2 = plot(lon1,lat1,'go');
    elseif strcmp(fininp,'e')
        if exist('pl3','var')
            delete(pl3);
            delete(pl4);
        end
        [x,y] = ginput(1);
        dd = sqrt((lon-x).^2+(lat-y).^2);
        [~,Imax] = min(dd);
        ddmax = sqrt((lat(Imax)-lat(idx)).^2+(lon(Imax)-lon(idx)).^2);
        lat2 = mean(lat(idx(ddmax<0.5)));
        lon2 = mean(lon(idx(ddmax<0.5)));
        pl3 = plot(lon(Imax),lat(Imax),'ro');
        pl4 = plot(lon2,lat2,'go');
    end
    if ~isnan(lat1) && ~isnan(lat2)
        if exist('pl5','var')
            delete(pl5);
        end
        az = azimuth(lat1, lon1, lat2, lon2);
        [latn, lonn, r] = gc2sc(lat1, lon1, az);
        [latc, lonc] = scircle1(latn, lonn, r,[],[],[],1000);
        ddmin = sqrt((lat(Imin)-latc).^2+(lon(Imin)-lonc).^2);
        ddmax = sqrt((lat(Imax)-latc).^2+(lon(Imax)-lonc).^2);
        [~,I1] = min(ddmin);
        [~,I2] = min(ddmax);
        if I1 > I2
            it = I1;
            I1 = I2;
            I2 = it;
        end
        kmlen = deg2km(distance(latc(I1-3),lonc(I1-3),latc(I2+3),lonc(I2+3)));
        sr = 0.1;
        [lato,lono] = gcwaypts(latc(I1-3),lonc(I1-3),latc(I2+3),lonc(I2+3),round(kmlen./sr));
        pl5 = plot(lono,lato,'r','LineWidth',1.1);
    end
end

foundflag = zeros(size(lato));
for i = 1:length(lato)
    inp(i).lat = lato(i);
    inp(i).lon = lono(i);
    inp(i).baz = 0;
    inp(i).dist = 0;
    inp(i).si = 0;
    inp(i).err = 0;
    inp(i).per = 0;
end

for i = 1:length(data)
    latd(i) = data(i).lat(1);
    lond(i) = data(i).lon(1);
end

n = 1;
for i = 1:length(data)
    dd = delaz(latd(i),lond(i),lato,lono,0);
    [~,I] = min(dd);
    if dd(I) > dmax
        disidx(n) = i; 
        n = n+1;
        continue
    end
    for j = 1:length(data(i).ev)
        if ~foundflag(I)
            inp(I).lato = data(i).lat(1);
            inp(I).lono = data(i).lon(1);
            inp(I).baz = data(i).ev(j).baz;
            inp(I).dist = data(i).ev(j).dist;
            inp(I).dep = data(i).ev(j).dep;
            inp(I).si = data(i).ev(j).si;
            inp(I).err = data(i).ev(j).sierr;
            inp(I).per = data(i).ev(j).per;
            foundflag(I) = 1;
        else
            inp(I).lato = [inp(I).lato data(i).lat(1)];
            inp(I).lono = [inp(I).lono data(i).lon(1)];
            inp(I).baz = [inp(I).baz data(i).ev(j).baz];
            inp(I).dist = [inp(I).dist data(i).ev(j).dist];
            inp(I).dep = [inp(I).dep data(i).ev(j).dep];
            inp(I).si = [inp(I).si data(i).ev(j).si];
            inp(I).err = [inp(I).err data(i).ev(j).sierr];
            inp(I).per = [inp(I).per data(i).ev(j).per];
        end
    end
% disp([num2str(i) ' of ' num2str(length(data))]);
end
inp(~foundflag) = [];
n = 1;
for i = 1:length(inp)
    ii = abs(inp(i).si)>3.5;
    inp(i).baz(ii) = [];
    inp(i).dist(ii) = [];
    inp(i).si(ii) = [];
    inp(i).err(ii) = [];
    inp(i).per(ii) = [];
    inp(i).dep(ii) = [];
    if isempty(inp(i).si)
        todel(n) = i;
        n = n+1;
    end
end
if exist('todel','var')
inp(todel) = [];
end

lat0 = inp(1).lat;
lon0 = inp(1).lon;
inp(1).x = 0;
xx(1) = 0;
for i = 2:length(inp)
    [dist,az(i)] = delaz(lat0,lon0,inp(i).lat,inp(i).lon,0);
    inp(i).x = deg2km(dist);
    xx(i) = inp(i).x;
end
caz = mean(az(2:end));
for i = 1:length(inp)
    inp(i).caz = caz;
end

for i = 1:length(inp)
    latp(i) = inp(i).lat;
    lonp(i) = inp(i).lon; 
    nn(i) = length(inp(i).si);  
end
fig = figure; scatter(lonp,latp,nn.*2,nn,'filled');
hold on;
plot(lond,latd,'k.')
if exist('disidx','var')
plot(lond(disidx),latd(disidx),'ro')
end
% plot(lonp,latp,'w*','LineWidth',1.2); 
% plot(lonp,latp,'k.');
kmlen = deg2km(distance(latp(1),lonp(1),latp(end),lonp(end)));
[latp2,lonp2] = gcwaypts(latp(1),lonp(1),latp(end),lonp(end),round(kmlen./sr));
plot(lonp2,latp2,'r','LineWidth',1.5)
indx = round(linspace(1,round(round(kmlen./sr)./500).*500,round(round(kmlen./sr)./500)+1));
indx2 = round(linspace(1,round(round(kmlen./sr)./1000).*1000,round(round(kmlen./sr)./1000)+1));
indx3 = round(linspace(1,round(round(kmlen./sr)./2000).*2000,round(round(kmlen./sr)./2000)+1));
indx4 = round(linspace(1,round(round(kmlen./sr)./5000).*5000,round(round(kmlen./sr)./5000)+1));
indx5 = round(linspace(1,round(round(kmlen./sr)./10000).*10000,round(round(kmlen./sr)./10000)+1));
indx6 = round(linspace(1,round(round(kmlen./sr)./25000).*25000,round(round(kmlen./sr)./25000)+1));
strpl = 50;
if length(indx) > 10
    indx = indx2;
    strpl = 100;
end
if length(indx) > 10
    indx = indx3;
    strpl = 200;
end
if length(indx) > 10
    indx = indx4;
    strpl = 500;
end
if length(indx) > 10
    indx = indx5;
    strpl = 1000;
end
if length(indx) > 10
    indx = indx6;
    strpl = 2500;
end
if indx(end) > length(lonp2)
    indx(end) = [];
end
plot(lonp2(indx),latp2(indx),'r.','MarkerSize',25); plot(lonp2(indx),latp2(indx),'w.','MarkerSize',18)

cb = colorbar; axis equal; axis([min(lon)-0.5 max(lon)+0.5 min(lat)-0.5 max(lat)+0.5]); hold on; 
xlabel('Longitude in [deg]')
ylabel('Latitude in [deg]')
ylabel(cb,'Number of events')
title(['Total number of events (' num2str(sum(nn)) ') ' num2str(strpl) ' km'])
if ~exist([sel_data.SI_dir '/graphics/'],'dir')
    mkdir([sel_data.SI_dir '/graphics/']);
end
grid on
print(fig,[sel_data.SI_dir '/graphics/StationMap.jpg'],'-r600','-djpeg')
close(fig)


for i = 1:length(inp)
    disp(['Preparing virtual station ' num2str(i) ' of ' num2str(length(inp))])
    for j = 1:length(inp(i).si)
        [phases]=get_q(inp(i).dist(j),inp(i).dep(j),psflag);
        if isstruct(phases)
            pp = phases(1).q;
        else
            pp = nan;
        end
%         disp(['station ' num2str(i) ' ' num2str(inp(i).dist(j)) ' ' num2str(pp)])
        inp(i).p(j) = pp;
    end
end
save([sel_data.SI_dir '/input/input2D.mat'],'inp');

end