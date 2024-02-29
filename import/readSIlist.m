function RES = readSIlist(sel_data)
% reading data from list

% Copyright 2024 F.Link and M.D.Long 

sel_data.SI_file
fileID = fopen(sel_data.SI_file,'rt');
C = textscan(fileID,'%s|%s|%f|%f|%f|%f|%f|%f|%f|%f','HeaderLines',1);
fclose(fileID);

stat = C{1};
net = C{1};
lat = C{1};
lon = C{1};
si = C{1};
sierr = C{1};
per = C{1};
dist = C{1};
dep = C{1};
baz = C{1};



n = 1;

% get station(s)
stations = unique(stat);   

% loop over stations
for mm = 1:length(stations)
    selected_station = char(stations{mm});
    st_ind = strcmp(char(selected_station),stat);

    disp(['Doing ' selected_station ' ...'])
        
% loop over all events
    for oo = 1:length(st_ind)
        if oo == 1
            RES(n).stat = char(stat{st_ind(oo)});
            RES(n).net = char(net{st_ind(oo)});
            RES(n).lat = lat(st_ind(oo));
            RES(n).lon = lon(st_ind(oo));
        end
% read splitting intensity
        RES(n).ev(m).baz = baz(st_ind(oo));
        RES(n).ev(m).dist = dist(st_ind(oo));
        RES(n).ev(m).dep = dep(st_ind(oo));
        RES(n).ev(m).per = per(st_ind(oo));
        RES(n).ev(m).si = si(st_ind(oo));
        RES(n).ev(m).sierr = sierr(st_ind(oo));
        m = m+1;
        
% END loop over events
    end

    n = n+1;
% END loop over stations
end

    
end