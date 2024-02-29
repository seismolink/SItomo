function prep3D(sel_data)
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
end
dmax = 0.6;


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
    fininp = input('Define upper left corner (u), lower left corner (l) or finish selection (y): ','s');
    if strcmp(fininp,'u')
        if exist('pl1','var')
            delete(pl1);
        end
        [lon1,lat1] = ginput(1);
        pl1 = plot(lon1,lat1,'go');
    elseif strcmp(fininp,'l')
        if exist('pl2','var')
            delete(pl2);
        end
        [lon2,lat2] = ginput(1);
        pl2 = plot(lon2,lat2,'ro');
    end
    if ~isnan(lat1) && ~isnan(lat2)
        if exist('pl3','var')
            delete(pl3);
        end
        pl3 = plot([lon1 lon2 lon2 lon1 lon1],[lat1 lat1 lat2 lat2 lat1],'r','LineWidth',1.1);
    end
end

for i = 1:length(data)
    latd(i) = data(i).lat(1);
    lond(i) = data(i).lon(1);
end

n = 1;
m = 1;
nn = 0;
for i = 1:length(data)
    inflag = latd(i) <= lat1 && latd(i) >= lat2 && lond(i) >= lon1 && lond(i) <= lon2;
    if ~inflag
        disidx(n) = i; 
        n = n+1;
        continue
    end
    inp(m).lato = data(i).lat(1);
    inp(m).lono = data(i).lon(1);
    o = 1;
    for j = 1:length(data(i).ev)
        if isempty(data(i).ev(j).si)
            continue
        end
        inp(m).baz(o) = data(i).ev(j).baz;
        inp(m).dist(o) = data(i).ev(j).dist;
        inp(m).dep(o) = data(i).ev(j).dep;
        inp(m).si(o) = data(i).ev(j).si;
        inp(m).err(o) = data(i).ev(j).sierr;
        inp(m).per(o) = data(i).ev(j).per;
        nn = nn+1;
        o = o+1;
    end
    m = m+1;
end

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

lat0 = lat2;
lon0 = mean([lon2 lon1]);
xx = zeros(size(inp));
yy = zeros(size(inp));
for i = 1:length(inp)
    [dist,~] = delaz(lat0,lon0,inp(i).lato,lon0,0);
    inp(i).y = deg2km(dist);
    yy(i) = inp(i).y;
    [dist,~] = delaz(inp(i).lato,lon0,inp(i).lato,inp(i).lono,0);
    if lon0<inp(i).lono
        inp(i).x = deg2km(dist);
    else
        inp(i).x = -deg2km(dist);
    end
    xx(i) = inp(i).x;
end
% caz = mean(az(2:end));
caz = 0;
for i = 1:length(inp)
    inp(i).caz = caz;
end

xx2 = zeros(size(data));
yy2 = zeros(size(data));
for i = 1:length(data)
    [dist,~] = delaz(lat0,lon0,data(i).lat,lon0,0);
    if data(i).lat>=lat0
        yy2(i) = deg2km(dist);
    else
        yy2(i) = -deg2km(dist);
    end
    [dist,~] = delaz(data(i).lat,lon0,data(i).lat,data(i).lon,0);
    if lon0<=data(i).lon
        xx2(i) = deg2km(dist);
    else
        xx2(i) = -deg2km(dist);
    end
end

xx3 = zeros(size(plon));
yy3 = zeros(size(plat));
for i = 1:length(plat)
    [dist,~] = delaz(lat0,lon0,plat(i),lon0,0);
    if plat(i)>=lat0
        yy3(i) = deg2km(dist);
    else
        yy3(i) = -deg2km(dist);
    end
end
for i = 1:length(plon)
    [dist,~] = delaz(lat0,lon0,lat0,plon(i),0);
    if lon0<plon(i)
        xx3(i) = deg2km(dist);
    else
        xx3(i) = -deg2km(dist);
    end
end

blat = [lat1 lat2];
blon = [lon1 lon2];
xx4 = zeros(size(blon));
yy4 = zeros(size(blat));
for i = 1:length(blat)
    [dist,~] = delaz(lat0,lon0,blat(i),lon0,0);
    if blat(i)>=lat0
        yy4(i) = deg2km(dist);
    else
        yy4(i) = -deg2km(dist);
    end
end
for i = 1:length(blon)
    [dist,~] = delaz(lat0,lon0,lat0,blon(i),0);
    if lon0<blon(i)
        xx4(i) = deg2km(dist);
    else
        xx4(i) = -deg2km(dist);
    end
end

% for i = 1:length(inp)
%     latp(i) = inp(i).lat;
%     lonp(i) = inp(i).lon; 
%     nn(i) = length(inp(i).si);  
% end
fig = figure; 
s = imagesc(xx3,yy3,rho);
set(gca,'Ydir','normal')
colormap(jet)
cb = colorbar;
ylabel(cb,'Number of Events')
hold on;
plot(xx2,yy2,'w.','MarkerSize',15)
plot(xx2,yy2,'k.')
if exist('disidx','var')
plot(xx2(disidx),yy2(disidx),'rx','LineWidth',1.1)
end
plot([min(xx4) max(xx4) max(xx4) min(xx4) min(xx4)],[min(yy4) min(yy4) max(yy4) max(yy4) min(yy4)],'r','LineWidth',1.1);
plot([mean(xx4) mean(xx4)],[min(yy4) max(yy4)],'w','LineWidth',2)
plot([min(xx4) max(xx4)],[0 0],'w','LineWidth',2)
plot([mean(xx4) mean(xx4)],[min(yy4) max(yy4)],'k','LineWidth',1.1)
plot([min(xx4) max(xx4)],[0 0],'k','LineWidth',1.1)

kmlen = max(xx4)-min(xx4);
dkmi = [50 100 200 500 1000 2500];
ii = 1;
Nkm = round(kmlen./dkmi(ii));
while Nkm > 10
    ii = ii+1;
    Nkm = round(kmlen./dkmi(ii));
end
strpl = num2str(dkmi(ii));
xx5 = 0:dkmi(ii):max(xx4);
xx5b = -xx5;
yy5 = 0:dkmi(ii):max(yy4);
yy5b = -yy5;
plot(zeros(size(yy5)),yy5,'k.','MarkerSize',25); plot(zeros(size(yy5)),yy5,'w.','MarkerSize',18)
plot(zeros(size(yy5b)),yy5b,'k.','MarkerSize',25); plot(zeros(size(yy5b)),yy5b,'w.','MarkerSize',18)
plot(xx5,min(yy4).*ones(size(xx5)),'k.','MarkerSize',25); plot(xx5,min(yy4).*ones(size(xx5)),'w.','MarkerSize',18)
plot(xx5b,min(yy4).*ones(size(xx5b)),'k.','MarkerSize',25); plot(xx5b,min(yy4).*ones(size(xx5b)),'w.','MarkerSize',18)

cb = colorbar; axis equal; axis([min(xx2)-0.2.*dkmi(ii) max(xx2)+0.2.*dkmi(ii) min(yy2)-0.2.*dkmi(ii) max(yy2)+0.2.*dkmi(ii)]); hold on; 
xlabel('x(EW) in [km]')
ylabel('y(NS) in [km]')
ylabel(cb,'Number of events')
title(['Total number of events (' num2str(sum(nn)) ') ' num2str(strpl) ' km'])
if ~exist([sel_data.SI_dir '/graphics/'],'dir')
    mkdir([sel_data.SI_dir '/graphics/']);
end
grid on

print(fig,[sel_data.SI_dir '/graphics/StationMap3D.jpg'],'-r600','-djpeg')
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
save([sel_data.SI_dir '/input/input3D.mat'],'inp');

end