function start_SItomo(varargin)
% Start function for Splitting Intensity Tomography

% optional input
% 'version':      Print the current version of the program
% 'check':        Check if there is a newer version available
% 'synthetics:    Produce synthetic data (experimental)
% 'import':       import SI-data
% 'prep2D':       prepare SI-data to constuct 2D Profile
% 'prep3D':       prepare SI-data as 3D array
% 'invert2Dgd':   invert Splitting Intensities with Gradient Descent in 2D
% 'invert2Dbfgs': invert Splitting Intensities with BFGS Algorithm in 2D
% 'invert2Dmcmc': invert Splitting Intensities with rj-McMC Algorithm in 2D
% 'invert3Dgd':   invert Splitting Intensities with Gradient Descent in 3D
% 'invert3Dbfgs': invert Splitting Intensities with BFGS Algorithm in 3D
% 'invert3Dmcmc': invert Splitting Intensities with rj-McMC Algorithm in 3D
% 'plot2D':       Process and plot results from 2D inversion
% 'plot3D':       Process and plot results from 2D inversion

% additional input
%   sel_data:          use predefined input to avoid invoking user input

% Copyright 2024 F.Link and M.D.Long 

close all
clc
% clear all
clearvars -except varargin
clear functions
clear global

curr_version = '1.0.0';
program_url = 'https://github.com/obspy/obspy/';

% add all sub directories to search path
addpath('./defaults/')
addpath('./graphics/')
addpath('./import/')
addpath('./inversion/')
addpath(pwd)

% Check if version control is requested or predefined input exists
optflag = 0;
preinput_flag = 0;
chck_flag = 0;
vrsn_flag = 0;
for i = 1:length(varargin)
    if strcmp(varargin(i),'version')
        vrsn_flag = 1;
        optflag = 1;
        disp(['Current version of SItomo: ' curr_version]);
        return
    end
    if strcmp(varargin(i),'check')
        chck_flag = 1;
        optflag = 1;
        disp(['Current version of SItomo: ' curr_version]);
        disp('Checking online for recent version. Please wait...')
        if ispc
            request_url = 'wget64';
        else
            request_url = 'wget';
        end
        request_url = strcat(request_url,' "',program_url,'releases/latest" -O temp.txt');
        [~,~] = system(request_url);
        fid = fopen('temp.txt','rt');
        line = fgetl(fid);
        while ischar(line)
            if contains(line,'<title>')
                disp(line(10:end-8));
                if ~contains(line,curr_version)
                    fprintf(strcat('It is likely that a new version is available.\n',...
                        'Visit %s to download the newest version.\n'),program_url);
                else
                    disp('Current version is the most recent.')
                end
                break;
            end
            line = fgetl(fid);
        end
        fclose(fid);
        delete('temp.txt')
        return
    end
    if strcmp(varargin(i),'remote')
        sel_data.SI_dir = char(varargin(i+1));
        sel_data.addname = char(varargin(i+2));
        sel_data.remoteflag = 1;
        preinput_flag = 1;
        break;
    end
end
if ~preinput_flag
    sel_data = 0;
end

rng('shuffle');

warning('off','all')

txt = strvcat({'SItomo - Copyright 2022 F.Link and M.D.Long',...
    'This program is free software: you can redistribute it and/or ',...
    'modify it under the terms of the GNU General Public License as ',...
    'published by the Free Software Foundation, either version 3 of ',...
    'the License, or (at your option) any later version.', ' ',...
    'This program is distributed in the hope that it will be useful, ',...
    'but WITHOUT ANY WARRANTY; without even the implied warranty of ',...
    'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the ',...
    'GNU General Public License for more details. ', ' ',...
    'You should have received a copy of the GNU General Public License ',...
    'along with this program. ',...
    'If not, see <http://www.gnu.org/licenses/>.'});
disp(txt);

% Check which actions are requested
syn_flag = 0;
import_flag = 0;
prep2D_flag = 0;
prep3D_flag = 0;
invrt2Dgd_flag = 0;
invrt2Dbfgs_flag = 0;
invrt2Dmcmc_flag = 0;
invrt3Dgd_flag = 0;
invrt3Dbfgs_flag = 0;
invrt3Dmcmc_flag = 0;
plt2D_flag = 0;
plt3D_flag = 0;
for i = 1:length(varargin)
    if strcmp(varargin(i),'synthetics')
        syn_flag = 1;
        optflag = 1;
    end
    if strcmp(varargin(i),'import')
        import_flag = 1;
        optflag = 1;
    end
    if strcmp(varargin(i),'prep2D')
        prep2D_flag = 1;
        optflag = 1;
    end
    if strcmp(varargin(i),'prep3D')
        prep3D_flag = 1;
        optflag = 1;
    end
    if strcmp(varargin(i),'invert2Dgd')
        invrt2Dgd_flag = 1;
        optflag = 1;
    end
    if strcmp(varargin(i),'invert2Dbfgs')
        invrt2Dbfgs_flag = 1;
        optflag = 1;
    end
    if strcmp(varargin(i),'invert2Dmcmc')
        invrt2Dmcmc_flag = 1;
        optflag = 1;
    end
    if strcmp(varargin(i),'invert3Dgd')
        invrt3Dgd_flag = 1;
        optflag = 1;
    end
    if strcmp(varargin(i),'invert3Dbfgs')
        invrt3Dbfgs_flag = 1;
        optflag = 1;
    end
    if strcmp(varargin(i),'invert3Dmcmc')
        invrt3Dmcmc_flag = 1;
        optflag = 1;
    end
    if strcmp(varargin(i),'plot2D')
        plt2D_flag = 1;
        optflag = 1;
    end
    if strcmp(varargin(i),'plot3D')
        plt3D_flag = 1;
        optflag = 1;
    end
end
if ~optflag
    disp(' ')
    disp(' ')
    disp('SItomo is a program for inverting Splitting Intensity measurements from')
    disp('teleseismic shear-wave splitting analysis for 2D and 3D structure.')
    disp('------------------------------------------------------------------------')
    disp('Basic Usage of SItomo: ')
    disp('start_SItomo([''option''])')
    disp(' ')
    disp('The following options are available:')
    disp('version        - Print the current version of the program.')
    disp('check          - Check if there is a newer version available.')
    disp('synthetics     - Produce synthetic data (experimental).')
    disp('import         - import SI-data.')
    disp('prep2D         - import SI-data as 2D array.')
    disp('prep3D         - import SI-data as 3D array.')
    disp('invert2Dgd     - invert Splitting Intensities with Gradient Descent in 2D.')
    disp('invert2Dbfgs   - invert Splitting Intensities with BFGS Algorithm in 2D.')
    disp('invert2Dmcmc   - invert Splitting Intensities with rj-McMC Algorithm in 2D.')
    disp('invert3Dgd     - invert Splitting Intensities with Gradient Descent in 3D.')
    disp('invert3Dbfgs   - invert Splitting Intensities with BFGS Algorithm in 3D.')
    disp('invert3Dmcmc   - invert Splitting Intensities with rj-McMC Algorithm in 3D.')
    disp('plot2D         - Process and plot results from 2D inversion.')
    disp('plot3D         - Process and plot results from 2D inversion.')
    disp('remote         - Run inversion with predefined settings.')
end

% perform requested actions in correct order
if syn_flag
    [sel_data] = makesyn(sel_data);
end
if import_flag
    [sel_data] = importSI(sel_data);
end
if prep2D_flag
    prep2D(sel_data);
end
if prep3D_flag
    prep3D(sel_data);
end
if invrt2Dgd_flag
    [sel_data] = invert2Dgd(sel_data);
end
if invrt2Dbfgs_flag
    [sel_data] = invert2Dbfgs(sel_data);
end
if invrt2Dmcmc_flag
    [sel_data] = invert2Dmcmc(sel_data);
end
if invrt3Dgd_flag
    [sel_data] = invert3Dgd(sel_data);
end
if invrt3Dbfgs_flag
    [sel_data] = invert3Dbfgs(sel_data);
end
if invrt3Dmcmc_flag
    [sel_data] = invert3Dmcmc(sel_data);
end
if plt2D_flag
    plotresults2D(sel_data);
end
if plt3D_flag
    plotresults3D(sel_data);
end


end