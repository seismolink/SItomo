function [sel_data] = import3D(sel_data)
% control function to gather nescessary parameters and invoke the
% import of SI data from SplitRacer for 3D inversion
%
% Copyright 2024 F.Link and M.D.Long  

if ~isstruct(sel_data)
    input_flag = 0;
    clear sel_data
else
    input_flag = 1;
end

if ~input_flag
    sel_data.SI_dir = input('Define working directory for SI tomography: ','s');
end

if ~input_flag
    finflag = 0;
    n = 1;
    while ~finflag
        sel_data.work_dir{n} = input('Define working directory from SplitRacer: ','s');
        addflag = input('Add another folder structure? 1=yes,0=no :    ');
        if ~addflag
            finflag = 1;
        end
    end
end



for n=1:length(sel_data.work_dir)
% Check if a settings file already exists and prepare parameters
A = dir([sel_data.work_dir{n},'/results/Filter*']);
if isempty(A)
    disp('No Single Splitting data available.')
    return
else
    disp(' ')
    disp('The following settings are available.')
    sfl = 1;
    for i = 1:length(A)
        Btemp = dir([sel_data.work_dir{n},'/results/',A(i).name,'/saved_settings.txt']);
        for j = 1:length(Btemp)
            B(sfl).folder = [sel_data.work_dir{n},'/results/',A(i).name,'/'];
            B(sfl).name = Btemp(j).name;
            disp([num2str(sfl) ' ' B(sfl).folder '/' B(sfl).name]);
            sfl = sfl+1;
        end
    end
    disp(' ')
    sfn = input('Choose a setting file for analysis: ');
    disp(['Chosen setting is: ' B(sfn).folder '/' B(sfn).name]);
    disp(' ')
    sel_data.results_folder{n} = B(sfn).folder;

    fileID = fopen([sel_data.results_folder{n},'saved_settings.txt']);
    C2 = textscan(fileID,'%f %f %3.1f %f %f %s');
    fclose(fileID);

    sel_data.p1(n) = C2{1};
    sel_data.p2(n) = C2{2};
    sel_data.snr(n) = C2{3};
    sel_data.NL(n) = C2{4};
    sel_data.fa_area(n) = C2{5};
    % call set_folder back to memory
    sel_data.set_folder{n} = [sel_data.work_dir{n},char(C2{6})];
end

if ~input_flag && ~isfield(sel_data,'cat_phases')
    sel_data.cat_phases = input('Use all (1) or only null and good measurements (0)? ');
end

import_SR(sel_data);
prep3D(sel_data);

end