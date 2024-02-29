function [phases]=get_q(eq_dist,eq_depth,psflag)

% table look up of ray parameters
% usage:
% eq_dist: distance earthquake - station
% eq_depth: event depth

% Copyright 2024 F.Link and M.D.Long 

% load ray parameter files
if psflag
load('PS.mat')
else
load('SKS.mat')
% load('SKKS.mat')
% load('SKIKS.mat')
% load('PKS.mat')
% load('PKIKS.mat')
end

% table look up with interpolation
if psflag
[q.PS]=interp_q(PS,eq_dist,eq_depth);
else
[q.SKS]=interp_q(SKS,eq_dist,eq_depth);
% [q.SKKS]=interp_q(SKKS,eq_dist,eq_depth);
% [q.SKIKS]=interp_q(SKIKS,eq_dist,eq_depth);
% [q.PKS]=interp_q(PKS,eq_dist,eq_depth);
% [q.PKIKS]=interp_q(PKIKS,eq_dist,eq_depth);
end

fn = fieldnames(q);

% remove non existent phases
for iF = 1:length(fn)
    af = fn{iF};
    if isnan(q.(af))
        q = rmfield(q,af);
    end
end


fn2 = fieldnames(q);
% rewrite fields
for iF2 = 1:length(fn2)
    af = fn2{iF2};
    phases(iF2).name = af;
    phases(iF2).q = q.(af);
end 

if exist('phases','var') == 0
    phases = 0;
end

end