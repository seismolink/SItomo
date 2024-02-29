function [model,result] = DoInvert3Dbfgs(data,thetaflag,vflag,wflag,zlim,nx,ny,nz,ix,iy,iz,edgeadd,sets,N1,N2,n1,n2,nevmax,initflag,strength_in,phi_in,theta_in,alpha,dec,showflag,printflag,graphicsfolder,filename,addname)
% Copyright 2024 F.Link and M.D.Long 

[workdir,~,~] = fileparts(graphicsfolder(1:end-1));
if ~showflag
    fid = fopen([workdir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Initialize results variables.\n');
    fclose(fid);
end
% initialize results variables
if ix~=nx || iy~=ny || iz~=nz
    aa = zeros(ix,iy,iz,N1);
    bb = aa;
    cc = aa;
else
    aa = zeros(nx,ny,nz,N1);
    bb = aa;
    cc = aa;
end
mf = zeros(N1,1);

if ~showflag
    fid = fopen([workdir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Prepare weights.\n');
    fclose(fid);
end
% Prepare weights
X_ST=data.x.*1000;
Y_ST=data.y.*1000;
x=linspace(min(X_ST),max(X_ST),nx+1);
y=linspace(min(Y_ST),max(Y_ST),ny+1);
ind = 1:length(X_ST);
if wflag
    w0 = zeros(size(data.x));
    for i = 1:nx
        for j = 1:ny
            stats(i).ind = ind(X_ST>=x(i)&X_ST<=x(i+1)&Y_ST>=y(j)&Y_ST<=y(j+1));
            w0(stats(i).ind) = 1./length(stats(i).ind);
        end
    end
else
    w0 = ones(size(data.x));
end

if ~showflag
    fid = fopen([workdir '/results/log' addname '.txt'],'at');
    fprintf(fid,'Start loop for statistics.\n');
    fclose(fid);
end

X_ST=data.x;
Y_ST=data.y;
if ix~=nx || iy~=ny || iz~=nz
    A_0 = zeros(ix,iy,iz);
    B_0 = zeros(ix,iy,iz);
    C_0 = zeros(ix,iy,iz);
else
    A_0 = zeros(nx,ny,nz);
    B_0 = zeros(nx,ny,nz);
    C_0 = zeros(nx,ny,nz);
end

% start loop for statistics
kk = 0;
while kk<N1
    % starting model
    if ix~=nx || iy~=ny || iz~=nz
        A_0 = zeros(ix,iy,iz);
        B_0 = zeros(ix,iy,iz);
        C_0 = zeros(ix,iy,iz);
    end
    if initflag
        A_0 = strength_in;
        B_0 = phi_in;
        C_0 = theta_in;
    else
        A_0(:) = 0.001;
        B_0(:) = rand.*pi;
        C_0(:) = 0;
    end
    
    if printflag
        if showflag
            fig = figure('Position',[300 50 600 580]);
        else
            fig = figure('Position',[300 50 600 580],'visible','off');
        end
    else
        fig = 0;
    end
    
    % use different subset for each statistical iteration
    noi = min(nevmax,round(length(data.x)*4/5));
    indoi = randperm(length(data.x),noi);
        
        if wflag
            w = zeros(size(data.x));
            for jj = 1:length(stats)
                stats(jj).ind2 = intersect(indoi,stats(jj).ind);
                w(stats(jj).ind2) = 1./length(stats(jj).ind2);
            end
        else
            w = ones(size(data.x));
        end
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Initialize input for matlab-internal function fminunc.\n');
            fclose(fid);
        end
        % calculate gradients
        funinp.data = data;
        funinp.indoi = indoi;
        funinp.weights = w;
        funinp.sflag = 0;
        funinp.vflag = vflag;
        funinp.zlim = zlim*1000;
        funinp.nx = nx;
        funinp.ny = ny;
        funinp.nz = nz;
        funinp.ix = ix;
        funinp.iy = iy;
        funinp.iz = iz;
        funinp.edge = edgeadd;
        funinp.sets = sets;
        funinp.thetaflag = thetaflag;
        funinp.showflag = showflag;
        funinp.initflag = initflag;
        funinp.dec = dec;
        funinp.filename = [workdir '/results/log' addname '.txt'];
        funinp.smooth = alpha;
        options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'MaxIterations',5,'MaxFunctionEvaluations',n2);
        problem.options = options;
        problem.x0 = [A_0(:)' B_0(:)' C_0(:)'];
        if initflag
            a0 = A_0(:);
            b0 = B_0(:);
            c0 = C_0(:);
            indx = a0~=0;
            problem.x0 = [a0(indx)' b0(indx)' c0(indx)'];
            funinp.indx = indx; 
        end
        problem.objective = @(x)SIsensitivityBFGS3D(x,funinp,fig);
        problem.solver = 'fminunc';
        for i = 1:n1
            if ~showflag
                fid = fopen([workdir '/results/log' addname '.txt'],'at');
                fprintf(fid,'Approximate starting model using short-cycles of BFGS-algorithm using matlab-internal function fminunc.\n');
                fclose(fid);
            end
            [mm,fvaltemp] = fminunc(problem);
            if ix~=nx || iy~=ny || iz~=nz
                if initflag
                    a0 = zeros(1,ix*iy*iz);
                    b0 = a0;
                    c0 = a0;
                    a0(indx) = mm(1:sum(indx));
                    b0(indx) = mm(sum(indx)+1:2*sum(indx));
                    c0(indx) = mm(2*sum(indx)+1:3*sum(indx));
%                     A_0 = reshape(a0(k)',nx,ny,nz);
%                     B_0 = reshape(b0(k)',nx,ny,nz);
%                     C_0 = reshape(c0(k)',nx,ny,nz);
                    A_0 = reshape(a0',ix,iy,iz);
                    B_0 = reshape(b0',ix,iy,iz);
                    C_0 = reshape(c0',ix,iy,iz);
                else
%                     A_0 = reshape(mm(k)',nx,ny,nz);
%                     B_0 = reshape(mm(k+max(k(:)))',nx,ny,nz);
%                     C_0 = reshape(mm(k+2.*max(k(:)))',nx,ny,nz);
                    A_0 = reshape(mm(1:ix*iy*iz)',ix,iy,iz);
                    B_0 = reshape(mm(ix*iy*iz+1:2*ix*iy*iz)',ix,iy,iz);
                    C_0 = reshape(mm(2*ix*iy*iz+1:3*ix*iy*iz)',ix,iy,iz);
                end
            else
                if initflag
                    a0 = zeros(1,nx*ny*nz);
                    b0 = a0;
                    c0 = a0;
                    a0(indx) = mm(1:sum(indx));
                    b0(indx) = mm(sum(indx)+1:2*sum(indx));
                    c0(indx) = mm(2*sum(indx)+1:3*sum(indx));
                    A_0 = reshape(a0,nx,ny,nz);
                    B_0 = reshape(b0,nx,ny,nz);
                    C_0 = reshape(c0,nx,ny,nz);
                else
                    A_0 = reshape(mm(1:nx*ny*nz),nx,ny,nz);
                    B_0 = reshape(mm(nx*ny*nz+1:2*nx*ny*nz),nx,ny,nz);
                    C_0 = reshape(mm(2*nx*ny*nz+1:3*nx*ny*nz),nx,ny,nz);
                end
            end
            ineg = A_0<0;
            A_0(ineg) = abs(A_0(ineg));
            B_0(ineg) = mod(B_0(ineg)+pi/2,pi);
            inull = A_0<0.005;
            A_0(inull) = 0.001;
            e1 = std(mod(B_0(inull),pi));
            temp = mod(B_0,pi);
            temp(temp>pi/2) = temp(temp>pi/2)-pi;
            e2 = std(temp(inull));
            if e1<e2
                B_0(inull) = mean(mod(B_0(inull),pi));
            else
                B_0(inull) = mean((temp(inull)));
            end
            C_0(inull) = 0;
%             if ix~=nx || iy~=ny || iz~=nz
%                 for j = 1:length(ku)
%                     A0n(j) = mean(A_0(k==ku(j)));
%                     B0n(j) = mean(B_0(k==ku(j)));
%                     C0n(j) = mean(C_0(k==ku(j)));
%                 end
%                 A_0 = reshape(A0n,ix,iy,iz);
%                 B_0 = reshape(B0n,ix,iy,iz);
%                 C_0 = reshape(C0n,ix,iy,iz);
%             end
            if initflag
                a0 = A_0(:);
                b0 = B_0(:);
                c0 = C_0(:);
                problem.x0 = [a0(indx)' b0(indx)' c0(indx)'];
            else
                problem.x0 = [A_0(:)' B_0(:)' C_0(:)'];
            end
        end
        if sqrt((fvaltemp./sum(abs(data.x(indoi))>0.2))*2)>0.5
            continue;
        else
            kk=kk+1;
        end
        options = optimoptions('fminunc','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'MaxIterations',N2);
        problem.options = options;
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Run inversion with BFGS-algorithm using matlab-internal function fminunc.\n');
            fclose(fid);
        end
        [mm,fval] = fminunc(problem);
        if ~showflag
            fid = fopen([workdir '/results/log' addname '.txt'],'at');
            fprintf(fid,'Step of the inversion finished.\n');
            fclose(fid);
        end
        
        
        if ix~=nx || iy~=ny || iz~=nz
            if initflag
                a0 = zeros(1,ix*iy*iz);
                b0 = a0;
                c0 = a0;
                a0(indx) = mm(1:sum(indx));
                b0(indx) = mm(sum(indx)+1:2*sum(indx));
                c0(indx) = mm(2*sum(indx)+1:3*sum(indx));
%                     A_0 = reshape(a0(k)',nx,ny,nz);
%                     B_0 = reshape(b0(k)',nx,ny,nz);
%                     C_0 = reshape(c0(k)',nx,ny,nz);
                A_0 = reshape(a0',ix,iy,iz);
                B_0 = reshape(b0',ix,iy,iz);
                C_0 = reshape(c0',ix,iy,iz);
            else
%                     A_0 = reshape(mm(k)',nx,ny,nz);
%                     B_0 = reshape(mm(k+max(k(:)))',nx,ny,nz);
%                     C_0 = reshape(mm(k+2.*max(k(:)))',nx,ny,nz);
                A_0 = reshape(mm(1:ix*iy*iz)',ix,iy,iz);
                B_0 = reshape(mm(ix*iy*iz+1:2*ix*iy*iz)',ix,iy,iz);
                C_0 = reshape(mm(2*ix*iy*iz+1:3*ix*iy*iz)',ix,iy,iz);
            end
        else
            if initflag
                a0 = zeros(1,nx*ny*nz);
                b0 = a0;
                c0 = a0;
                a0(indx) = mm(1:sum(indx));
                b0(indx) = mm(sum(indx)+1:2*sum(indx));
                c0(indx) = mm(2*sum(indx)+1:3*sum(indx));
                A_0 = reshape(a0,nx,ny,nz);
                B_0 = reshape(b0,nx,ny,nz);
                C_0 = reshape(c0,nx,ny,nz);
            else
                A_0 = reshape(mm(1:nx*ny*nz),nx,ny,nz);
                B_0 = reshape(mm(nx*ny*nz+1:2*nx*ny*nz),nx,ny,nz);
                C_0 = reshape(mm(2*nx*ny*nz+1:3*nx*ny*nz),nx,ny,nz);
            end
        end
        ineg = A_0<0;
        A_0(ineg) = abs(A_0(ineg));
        B_0(ineg) = mod(B_0(ineg)+pi/2,pi);
        inull = A_0<0.005;
        A_0(inull) = 0;
        e1 = std(mod(B_0(:),pi));
        temp = mod(B_0,pi);
        temp(temp>pi/2) = temp(temp>pi/2)-pi;
        e2 = std(temp(:));
        if e1<e2
            B_0(inull) = mean(mod(B_0(:),pi));
        else
            B_0(inull) = mean((temp(:)));
        end
        C_0(inull) = 0;
        
    eflag(:,kk) = zeros(size(data.x));
    eflag(indoi,kk) = 1;
    aa(:,:,:,kk) = A_0;
    bb(:,:,:,kk) = mod(B_0,pi);
    cc(:,:,:,kk) = C_0;
    mf(kk) = fval;
    if ~showflag
        fid = fopen([workdir '/results/log' addname '.txt'],'at');
        fprintf(fid,'Current step stored in variables.\n');
        fclose(fid);
    end
    if printflag
        filename2 = [filename '_' addname num2str(kk) '_.jpg'];
        print(fig,[graphicsfolder '/' filename2],'-djpeg','-r300');
        close(fig)
        fid = fopen([workdir '/results/log' addname '.txt'],'at');
        fprintf(fid,'Current plot printed to file.\n');
        fclose(fid);
    end
end

model.sets = sets;
model.zlim = zlim;
model.edges = edgeadd;
model.nx = nx;
model.ny = ny;
model.nz = nz;
model.ix = ix;
model.iy = iy;
model.iz = iz;
model.N1 = N1;
model.N2 = N2;
model.n1 = n1;
model.n2 = n2;
model.vflag = vflag;
model.thetaflag = thetaflag;
model.initflag = initflag;
model.initx = strength_in;
model.initphi = phi_in;
model.inittheta = theta_in;
model.data = data;
model.w = w0;
result.eflag = eflag;
result.x = aa;
result.phi = bb;
result.theta = cc;
result.mf = mf;

end