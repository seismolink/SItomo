function [delta,azeqst,azsteq]=delaz(eqlat,eqlon,stlat,stlon,flag)

% compute earthquake/station distance and azimuth

% usage: [delta,azeqst,azsteq]=delaz(eqlat,eqlon,stlat,stlon,flag);
%
%     compute distance and azimuth from earthquake (eq) to station (st)
%     delta  = distance between (eq) and (st) in degrees
%     azeqst = azimuth from (eq) to (st) clockwise from north in degrees
%     azsteq = azimuth from (st) to (eq) clockwise from north in degrees
%
%     if input coordinates are geographic degrees   flag=0
%     if input coordinates are geocentric radians   flag=1
%
%     input latitudes and longitudes can be scalars or vectors
%     acceptable combinations are one earthquake with many stations,
%     one station and many earthquakes, or the same number of earthquakes
%     and stations
%     output vectors will have same dimensions as input vectors
%
%     calls coortr.m

% convert from geographic degrees to geocentric radians if necessary
% convert to spherical polar coordinates in radians (lat -> colatitude)

% taken from the CORAL package of K. Creager, 
% http://earthweb.ess.washington.edu/creager/coral_tools.html
% ftp://ftp.geophys.washington.edu/pub/out/corat.tar.Z, 
% Seismological Research Letters, Volume 68, Number 2, 1997

if flag==0   % convert geographic degrees to geocentric radians
  [eqlat,eqlon]=coortr(eqlat,eqlon,flag); 
  [stlat,stlon]=coortr(stlat,stlon,flag); 
end

eqcolat=pi/2-eqlat;
stcolat=pi/2-stlat;

cos_eq=cos(eqcolat);
sin_eq=sin(eqcolat);
cos_st=cos(stcolat);
sin_st=sin(stcolat);
cos_eqst=cos(stlon-eqlon);
sin_eqst=sin(stlon-eqlon);

cos_delta=cos_eq.*cos_st + sin_eq.*sin_st.*cos_eqst;
sin_delta=sqrt(1-cos_delta.*cos_delta);
delta=atan2(sin_delta,cos_delta);

sin_delta = sin_delta + (sin_delta==0)*eps;

% index is zero if expression is false, 1 if true; 
% if false, leave unchanged, if true azeqst=pi-azeqst
% this puts azeqst into the correct quadrant
azeqst=asin(sin_st.*sin_eqst./sin_delta);
index=(sin_eq.*cos_st - cos_eq.*sin_st.*cos_eqst < 0);
azeqst=azeqst + index.*(pi-2*azeqst);
azeqst=azeqst + (azeqst<0)*2*pi;

azsteq=asin(-sin_eq.*sin_eqst./sin_delta);
index=(cos_eq.*sin_st - sin_eq.*cos_st.*cos_eqst < 0);
azsteq=azsteq + index.*(pi-2*azsteq);
azsteq=azsteq + (azsteq<0)*2*pi;

% convert to degrees
delta=delta*180/pi;
azeqst=azeqst*180/pi;
azsteq=azsteq*180/pi;
