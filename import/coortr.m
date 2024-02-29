function [latout,lonout]=coortr(latin,lonin,flag)

% geocentric/geographic coordinate transformation

% usage: [latout,lonout]=coortr(latin,lonin,flag);

%  purpose: transform between geographic and geocentric coordinates
%           geographic degrees to geocentric radians ( if flag=0 )
%        or geocentric radians to geographic degrees ( if flag=1 )
%           earthquake and station locations are typically
%           given in geographic coordinates, whereas most calculations
%           such as epicentral distance are given in geocentric coordinates
%           latin and lonin can be scalars or vectors,
%           latout and lonout will match the dimensions of latin and lonin

%  if flag==0,
%     latin,  lonin  are latitude and longitude in geographic degrees
%     latout, lonout are latitude and longitude in geocentric radians
%  if flag==1,
%     latin,  lonin  are latitude and longitude in geocentric radians
%     latout, lonout are latitude and longitude in geographic degrees


% taken from the CORAL package of K. Creager, 
% http://earthweb.ess.washington.edu/creager/coral_tools.html
% ftp://ftp.geophys.washington.edu/pub/out/corat.tar.Z, 
% Seismological Research Letters, Volume 68, Number 2, 1997

if flag==0
    latout=atan(tan(latin*pi/180)*0.9933056);
    lonout=lonin*pi/180;
elseif flag==1
    latout=atan(tan(latin)/0.9933056)*180/pi;
    lonout=lonin*180/pi;
end

end
