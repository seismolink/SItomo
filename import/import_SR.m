function RES = import_SR(sel_datatmp)
% reading data from SplitRacer structure

% Copyright 2024 F.Link and M.D.Long 

n = 1;

% loop over input folders
for nn = 1:length(sel_datatmp.results_folder)
    sel_data.station = 'all';
    sel_data.work_dir = sel_datatmp.work_dir{nn};
    sel_data.results_dir = sel_datatmp.results_folder{nn};
    sel_data.reevaluate = sel_datatmp.reevaluate;

% get station(s)
    [~,stations] = read_stationfile(sel_data);    

% loop over stations
    for mm = 1:length(stations.name)
        selected_station = char(stations.name{mm});
        st_ind = strcmp(char(selected_station),stations.name);
        dir_get = [sel_data.results_dir '/' char(stations.name{st_ind}) '_' char(stations.nc{st_ind}) '/'];
        
        disp(['Doing ' selected_station ' ...'])
        
% check if single splitting results exist
        st_get = char(strcat(dir_get,'/results_si_split.mat'));
        if ~exist(st_get,'file')
            continue
        end

% load single splitting results exist
        load(st_get);
        fn = fieldnames(results);
        m = 1;
        
% load categories
        filename = strcat(dir_get,'/phases_js.mat');
        if exist(filename,'file')
            load(filename);
        else
            continue
        end
        
% loop over all events
        for oo = 1:length(fn)
            if contains(fn{oo},'event')
                event = results.(fn{oo});
            else
                continue
            end
            
            foundflag = 0;
            fn2 = fieldnames(data);
            for j = 1:length(fn2)
                if ~isempty(strfind(fn2{j},'event'))
                    if data.(fn2{j}).origin_time==event.origin_time && data.(fn2{j}).i_phase == event.i_phase
                        cat = data.(fn2{j}).cat;
                        foundflag = 1;
                        break;
                    end
                end
            end
            if ~foundflag
                cat = 'poor';
            end  
            
            if m == 1
                RES(n).stat = char(stations.name{st_ind});
                RES(n).net = char(stations.nc{st_ind});
                RES(n).lat = stations.lat(st_ind);
                RES(n).lon = stations.lon(st_ind);
            end
% read splitting intensity
            RES(n).ev(m).baz = event.baz;
            RES(n).ev(m).dist = event.dist;
            RES(n).ev(m).dep = event.depth;
            RES(n).ev(m).cat = cat;
            if ~sel_data.reevaluate
                RES(n).ev(m).per = 10;
                RES(n).ev(m).si = event.split_int;
                RES(n).ev(m).sierr = diff(event.split_err);
            else
                tw1 = event.tw1;
                tw2 = event.tw2;
                rad = event.rad;
                tra = event.tra;
                trtime = linspace(tw1,tw2,length(rad));
                dtn2 = trtime(2)-trtime(1);
                disp('Starting analysis')
                [~,~,~,~,~,~,~,~,si1,p,~,~] = SplitIntFunction_fin(tra,rad,dtn2);
                disp('Analysis done')
                
                todel = isinf(p);
                si1(todel) = [];
                p(todel) = [];
                
                siest = mean(si1);
                pest = mean(p);
                sierr = std(si1);
                
                if siest>2.2
                    continue
                end
                
                RES(n).ev(m).si = siest;
                RES(n).ev(m).sierr = sierr;
                RES(n).ev(m).per = pest;
            end
            m = m+1;
        
% END loop over events
        end
        
        n = n+1;
% END loop over stations
    end
% END loop over input folders
end
    
end