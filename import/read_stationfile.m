function [selected_stations,stations]=read_stationfile(sel_data)

% function finds selected stations from pull-down menu
% usage: all necessary input is stored in sel_data

% Copyright 2016 M.Reiss and G.Rümpker

%read stations.lst file
fileID1 = fopen([sel_data.work_dir,'/station.lst']);
C = textscan(fileID1,'%s %s %f %f');
fclose(fileID1);

stations.name = C{1,1};
stations.nc = C{1,2};
stations.lat = C{1,3};
stations.lon = C{1,4};

%define which stations should be processed

find_stations = strcmp('all',char(sel_data.station));

if find_stations == 0
    selected_stations = sel_data.station;
elseif find_stations == 1
    selected_stations = stations.name;
end


end