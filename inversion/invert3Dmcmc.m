function [sel_data] = invert3Dmcmc(sel_data)
% inversion for 3D anisotropic structure based on splitting intensities
%
% Copyright 2024 F.Link and M.D.Long 

% This is the settings and call function for the inversion based on a 
% gradient informed reversible-jump Markov chain Monte Carlo with 
% Metropolis-Hastings algorithm.

if ~isstruct(sel_data) 
    input_flag = 0;
    clear sel_data
    sel_data.remoteflag = 0;
    addname = input('Define string ID for current inversion: ','s');
else
    addname = sel_data.addname;
    if ~sel_data.remoteflag
        input_flag = 0;
    else
        input_flag = 1;
    end
end

if ~input_flag
    sel_data.SI_dir = input('Define working directory for SI tomography: ','s');
end
if ~exist([sel_data.SI_dir '/input/input3D.mat'],'file')
    disp('No input structure. Import splitting intensity data first for 2D structure.')
    if sel_data.remoteflag
        fid = fopen([sel_data.SI_dir '/results/log' addname '.txt'],'at');
        fprintf(fid,'No input structure. Import splitting intensity data first for 2D structure.\n');
        fclose(fid);
    end
    return
else
    if sel_data.remoteflag
        fid = fopen([sel_data.SI_dir '/results/log' addname '.txt'],'at');
        fprintf(fid,'Load 3D input structure.\n');
        fclose(fid);
    end
    load([sel_data.SI_dir '/input/input3D.mat']);
end
if ~exist([sel_data.SI_dir '/results/'],'dir')
    mkdir([sel_data.SI_dir '/results/']);
end
if sel_data.remoteflag
    fid = fopen([sel_data.SI_dir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Gradient informed rj-McMC algorithm chosen for inversion.\n');
    fclose(fid);
end

if input_flag
    newflag = 0;
else
    newflag = 1;
    if exist([sel_data.SI_dir '/results/settings3D.txt'],'file')
        newflag = input('Settings file exists. Use existing (0) or define new settings (1)? :  ');
    end
end
if newflag
thetaflag = input('Search for dipping symmetry axis? 1=yes,0=no :  ');
vflag = input('Consider non-vertical incidence? 1=yes,0=no :   ');
if vflag == 1
    vflag = 0;
else
    vflag = 1;
end
z1 = input('Define minimum depth in [km]:  ');
z2 = input('Define maximum depth in [km]:  ');
nx = input('Define model dimension in x:  ');
ny = input('Define model dimension in y:  ');
nz = input('Define model dimension in z:  ');
edgeadd = input('Define width of boundary region in [km] (e.g., 100):  ');
edgeadd = edgeadd.*1000;
xstepini = input('Define initial step size for anisotropy (e.g., 0.25):  ');
phistepini = input('Define initial step size for fast axis in [deg] (e.g., 30):  ');
phistepini = phistepini/180*pi;
if ~vflag
    thetastepini = input('Define initial step size for dip in [deg] (e.g., 5):  ');
    thetastepini = thetastepini/180*pi;
else
    thetastepini = 0;
end
N1 = input('Define number of initial voronoi cells:  ');
nevmax = input('Define maximum number of events in the inversion: ');
N2 = input('Define number of iterations in inversion:  ');
n1 = input('Define after how many steps the reference should be re-calculated:  ');
n2 = input('Define minimum number of voronoi cells:  ');
initflag = input('Use initial model (requires init.mat in input folder)? 1=yes,0=no :  ');
showflag = input('Show intermediate results? 1=yes,0=no :  ');
printflag = input('Print intermediate results? 1=yes,0=no :  ');
maxfrac = 1;
alpha = input('Define damping parameter: ');
dec = input('Define distance weighting exponent: ');

% save settings
fid = fopen([sel_data.SI_dir '/results/settings3D.txt'],'wt');
fprintf(fid,'thetaflag|%f\n',thetaflag); % Theta as free parameter
fprintf(fid,'vflag|%f\n',vflag); % restricted to vertical incidence
fprintf(fid,'z1|%f\n',z1); % Minimum depth of model
fprintf(fid,'z2|%f\n',z2); % Maximum depth of model
fprintf(fid,'nx|%f\n',nx); % Number of Elements in x (underlying grid)
fprintf(fid,'nz|%f\n',nz); % Number of Elements in z (underlying grid)
fprintf(fid,'ny|%f\n',ny); % Number of Elements in y (>1 only in 3D)
fprintf(fid,'edgeadd|%f\n',edgeadd); % added volume at edges in meters
fprintf(fid,'xstepini|%f\n',xstepini); % Initial stepsize in x (only dd-algorithm)
fprintf(fid,'phistepini|%f\n',phistepini); % Initial stepsize in phi (only dd-algorithm)
fprintf(fid,'thetastepini|%f\n',thetastepini); % Initial stepsize in theta (only dd-algorithm)
fprintf(fid,'N1|%f\n',N1); % Number of initial Voronoi Cells
fprintf(fid,'nevmax|%f\n',nevmax); % Maximum number of events considered in inversion
fprintf(fid,'N2|%f\n',N2); % Number of Iterations in the McMC
fprintf(fid,'n1|%f\n',n1); % Re-calculate reference after n1 steps
fprintf(fid,'n2|%f\n',n2); % Minimum number of Voronoi Cells
fprintf(fid,'initflag|%f\n',initflag); % Use pre-defined starting model
fprintf(fid,'showflag|%f\n',showflag); % Prompt figure of intermediate steps
fprintf(fid,'printflag|%f\n',printflag); % Print figure at final step
fprintf(fid,'maxfrac|%f\n',maxfrac);
fprintf(fid,'alpha|%f\n',alpha); % Smoothing parameter (smooth > alpha > heterogeneous)
fprintf(fid,'dec|%f\n',dec); % Distance weighting from station (s ~ R.^dec)
fclose(fid);
else
    if sel_data.remoteflag
        fid = fopen([sel_data.SI_dir '/results/log' addname '.txt'],'at');
        fprintf(fid,'Read data for inversion from settings3D.txt\n');
        fclose(fid);
    end
    fid = fopen([sel_data.SI_dir '/results/settings3D.txt'],'rt');
    C = textscan(fid,'%s%f','Delimiter','|');
    fclose(fid);
    thetaflag = C{2}(1);
    vflag = C{2}(2);
    z1 = C{2}(3);
    z2 = C{2}(4);
    nx = C{2}(5);
    nz = C{2}(6);
    ny = C{2}(7);
    edgeadd = C{2}(8);
    xstepini = C{2}(9);
    phistepini = C{2}(10);
    thetastepini = C{2}(11);
    N1 = C{2}(12);
    nevmax = C{2}(13);
    N2 = C{2}(14);
    n1 = C{2}(15);
    n2 = C{2}(16);
    initflag = C{2}(17);
    showflag = C{2}(18);
    printflag = C{2}(19);
    maxfrac = C{2}(20);
    alpha = C{2}(21);
    dec = C{2}(22);
end

if initflag
    load([sel_data.SI_dir '/input/init.mat']);
    strength_in = finmod.x;
    phi_in = finmod.phiorig;
    theta_in = finmod.theta;
else
    strength_in = 0;
    phi_in = 0;
    theta_in = 0;
end

breakflag = 0;
n = 1;
for i = 1:length(inp)
    for j = 1:length(inp(i).si)
        if abs(inp(i).si(j))<2.0
            data.x(n) = inp(i).x;
            data.y(n) = inp(i).y;
            data.baz(n) = inp(i).baz(j);
            data.az(n) = 90-inp(i).baz(j);
            data.per(n) = inp(i).per(j);
            data.si(n) = inp(i).si(j);
            data.err(n) = inp(i).err(j);
            data.p(n) = inp(i).p(j);
            n = n+1;
        end
    end
    if breakflag
        break
    end
end
todela = isnan(data.si);
todelb = data.per>20;
todelc = data.per<4;
todel = (todela+todelb+todelc)>0;
data.x(todel) = [];
data.y(todel) = [];
data.baz(todel) = [];
data.az(todel) = [];
data.per(todel) = [];
data.p(todel) = [];
data.si(todel) = [];
data.err(todel) = [];

% calculate sensitivity
if ~exist([sel_data.SI_dir '/graphics/'],'dir')
    mkdir([sel_data.SI_dir '/graphics/']);
end
if sel_data.remoteflag
    fid = fopen([sel_data.SI_dir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Start inversion...\n');
    fclose(fid);
    starttime = tic;
end
[model,result] = DoInvert3Dmcmc(data,thetaflag,vflag,[z1 z2],nx,ny,nz,edgeadd,xstepini,phistepini,thetastepini,N1,N2,n1,n2,nevmax,initflag,strength_in,phi_in,theta_in,alpha,showflag,printflag,[sel_data.SI_dir '/graphics/'],'invert2Dres',addname);

save([sel_data.SI_dir '/results/result3D' addname '.mat'],'model','result');

if sel_data.remoteflag
    fid = fopen([sel_data.SI_dir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Inversion finished and succefully stored to result.mat.\n');
    fprintf(fid,'Inversion finished after %f seconds.\n', toc(starttime));
    fclose(fid);
end
end