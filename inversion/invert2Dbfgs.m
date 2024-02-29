function [sel_data] = invert2Dbfgs(sel_data)
% inversion for 2D anisotropic structure based on splitting intensities
%
% Copyright 2024 F.Link and M.D.Long 
%
% This is the settings and call function for the inversion based on an 
% Broyden–Fletcher–Goldfarb–Shanno (BFGS) Algorithm to estimate the Hessian
% matrix.

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
if ~exist([sel_data.SI_dir '/input/input2D.mat'],'file')
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
        fprintf(fid,'Load 2D input structure.\n');
        fclose(fid);
    end
    load([sel_data.SI_dir '/input/input2D.mat']);
end
if ~exist([sel_data.SI_dir '/results/'],'dir')
    mkdir([sel_data.SI_dir '/results/']);
end
if sel_data.remoteflag
    fid = fopen([sel_data.SI_dir '/results/log' addname '.txt'],'at');
    fprintf(fid,'BFGS-algorithm chosen for inversion.\n');
    fclose(fid);
end

if input_flag
    newflag = 0;
else
    newflag = 1;
    if exist([sel_data.SI_dir '/results/settings2D.txt'],'file')
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
wflag = input('Downweight stations with large amount of data? 1=yes,0=no :  ');
z1 = input('Define minimum depth in [km]:  ');
z2 = input('Define maximum depth in [km]:  ');
nx = input('Define model dimension in x:  ');
nz = input('Define model dimension in z:  ');
ix = input('Define free parameter dimension in x:  ');
iz = input('Define free parameter dimension in z:  ');
ny = 1;
edgeadd = input('Define width of boundary region in [km] (e.g., 100):  ');
edgeadd = edgeadd.*1000;
nny = 3;
N1 = input('Define number of statistical iterations:  ');
nevmax = input('Define maximum number of events per statistical iteration: ');
N2 = input('Define number of iterations in inversion:  ');
n1 = input('Define number of pre-cycles to prepare initial guess (e.g., 3):  ');
n2 = input('Define maximum number of iterations in pre-cycles (e.g., 10):  ');
initflag = input('Use initial model (requires init.mat in input folder)? 1=yes,0=no :  ');
showflag = input('Show intermediate results? 1=yes,0=no :  ');
printflag = input('Print intermediate results? 1=yes,0=no :  ');
maxfrac = 1;
alpha = input('Define damping parameter: ');
dec = input('Define magnitude of proximity penalty (exmpl: 0.5): ');

% save settings
fid = fopen([sel_data.SI_dir '/results/settings2D.txt'],'wt');
fprintf(fid,'thetaflag|%f\n',thetaflag); % Theta as free parameter
fprintf(fid,'vflag|%f\n',vflag); % restricted to vertical incidence
fprintf(fid,'wflag|%f\n',wflag); % weighting for amount of events per location
fprintf(fid,'z1|%f\n',z1); % Minimum depth of model
fprintf(fid,'z2|%f\n',z2); % Maximum depth of model
fprintf(fid,'nx|%f\n',nx); % Number of Elements in x
fprintf(fid,'nz|%f\n',nz); % Number of Elements in z
fprintf(fid,'ny|%f\n',ny); % Number of Elements in y (>1 only in 3D)
fprintf(fid,'ix|%f\n',ix); % Number of free parameters in x
fprintf(fid,'iz|%f\n',iz); % Number of free parameters in z
fprintf(fid,'edgeadd|%f\n',edgeadd); % added volume at edges in meters
fprintf(fid,'nny|%f\n',nny); % Resolution of y relative to dx (integer | 1 - dy=dx)
fprintf(fid,'N1|%f\n',N1); % Number of initial guesses (starting a new inversion)
fprintf(fid,'nevmax|%f\n',nevmax); % Maximum number of events considered in inversion
fprintf(fid,'N2|%f\n',N2); % Maximum number of Iterations
fprintf(fid,'n1|%f\n',n1); % Number of pre-cycles to prepare initial guess
fprintf(fid,'n2|%f\n',n2); % Maximum number of Function Evaluation in Pre-Cycles
fprintf(fid,'initflag|%f\n',initflag); % Use pre-defined starting model
fprintf(fid,'showflag|%f\n',showflag); % Prompt figure of intermediate steps
fprintf(fid,'printflag|%f\n',printflag); % Print figure at final step
fprintf(fid,'maxfrac|%f\n',maxfrac);
fprintf(fid,'alpha|%f\n',alpha); % Smoothing parameter (smooth > alpha > heterogeneous)
fprintf(fid,'dec|%f\n',dec); % Penalty weight for distance/depth (R^dec)
fclose(fid);
else
    if sel_data.remoteflag
        fid = fopen([sel_data.SI_dir '/results/log' addname '.txt'],'at');
        fprintf(fid,'Read data for inversion from settings2D.txt\n');
        fclose(fid);
    end
    fid = fopen([sel_data.SI_dir '/results/settings2D.txt'],'rt');
    C = textscan(fid,'%s%f','Delimiter','|');
    fclose(fid);
    thetaflag = C{2}(1);
    vflag = C{2}(2);
    wflag = C{2}(3);
    z1 = C{2}(4);
    z2 = C{2}(5);
    nx = C{2}(6);
    nz = C{2}(7);
    ny = C{2}(8);
    ix = C{2}(9);
    iz = C{2}(10);
    edgeadd = C{2}(11);
    nny = C{2}(12);
    N1 = C{2}(13);
    nevmax = C{2}(14);
    N2 = C{2}(15);
    n1 = C{2}(16);
    n2 = C{2}(17);
    initflag = C{2}(18);
    showflag = C{2}(19);
    printflag = C{2}(20);
    maxfrac = C{2}(21);
    alpha = C{2}(22);
    dec = C{2}(23);
end

if initflag
    load([sel_data.SI_dir '/input/init.mat']);
    strength_in = finmod.x;
    phi_in = finmod.phiorig;
    theta_in = finmod.theta;
    if nx ~= ix || nz ~= iz
        strength_in = reshape(strength_in,nx,nz);
        phi_in = reshape(phi_in,nx,nz);
        theta_in = reshape(theta_in,nx,nz);
        x = 1:nx; z = 1:nz;
        xn=linspace(1,nx,ix+1);
        zn=linspace(1,nz,iz+1);
        dxn = xn(2)-xn(1);
        dzn = zn(2)-zn(1);
        xn(end) = [];
        zn(end) = [];
        xn = xn+dxn/2;
        zn = zn+dzn/2;
        [X,Z] = meshgrid(x,z);
        [Xn,Zn] = meshgrid(xn,zn);
        Xn = Xn'; Zn = Zn'; X = X'; Z = Z';
        k = dsearchn([Xn(:) Zn(:)],[X(:) Z(:)]);
        ku = unique(k);
        for j = 1:length(ku)
            A0n(j) = mean(strength_in(k==ku(j)));
            B0n(j) = mean(phi_in(k==ku(j)));
            C0n(j) = mean(theta_in(k==ku(j)));
        end
        strength_in = reshape(A0n,ix,iz);
        phi_in = reshape(B0n,ix,iz);
        theta_in = reshape(C0n,ix,iz);
    end
else
    strength_in = 0;
    phi_in = 0;
    theta_in = 0;
end

n = 1;
for i = 1:length(inp)
    for j = 1:length(inp(i).si)
        data.x(n) = inp(i).x;
        data.y(n) = 0;
        data.baz(n) = inp(i).baz(j)+(90-inp(i).caz);
        data.az(n) = 90-inp(i).baz(j);
        data.per(n) = inp(i).per(j);
        data.si(n) = inp(i).si(j);
        data.err(n) = inp(i).err(j);
        data.p(n) = inp(i).p(j);
        n = n+1;
    end
end
todela = isnan(data.si);
todelb = data.per>20;
todelc = data.per<4;
todeld = isnan(data.per);
todel = (todela+todelb+todelc+todeld)>0;
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
if exist([sel_data.SI_dir '/input/sets2D.mat'],'file')
    load([sel_data.SI_dir '/input/sets2D.mat']);
else
    if sel_data.remoteflag
        fid = fopen([sel_data.SI_dir '/results/log' addname '.txt'],'at');
        fprintf(fid,'Get sampling from event distribution...\n');
        fclose(fid);
    end
    sets = getsampling(data,maxfrac,vflag,dec,[z1 z2],nx,ny,nz,edgeadd,nny,showflag,printflag,[sel_data.SI_dir '/graphics/'],'test');
    save([sel_data.SI_dir '/input/sets2D.mat'],'sets');
    if sel_data.remoteflag
        fid = fopen([sel_data.SI_dir '/results/log' addname '.txt'],'at');
        fprintf(fid,'Sampling file sets2D.mat created.\n');
        fclose(fid);
    end
end
if sel_data.remoteflag
    fid = fopen([sel_data.SI_dir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Start inversion...\n');
    fclose(fid);
end
[model,result] = DoInvert2Dbfgs(data,thetaflag,vflag,wflag,[z1 z2],nx,ny,nz,ix,iz,edgeadd,nny,sets,xstepini,phistepini,thetastepini,N1,N2,n1,n2,nevmax,initflag,strength_in,phi_in,theta_in,alpha,dec,showflag,printflag,[sel_data.SI_dir '/graphics/'],'invert2Dres',addname);

save([sel_data.SI_dir '/results/result' addname '.mat'],'model','result');

if sel_data.remoteflag
    fid = fopen([sel_data.SI_dir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Inversion finished and succefully stored to result.mat.\n');
    fclose(fid);
end
end