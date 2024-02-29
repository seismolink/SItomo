function [q_final]=interp_q(phase,eq_dist,eq_depth)

% table look up for travel times
% usage: 
% phase: phase name
% eq_dist: distance between station and event
% eq_depth: event depth

% Copyright 2024 F.Link and M.D.Long 

% define parameters - same matrix size as stored phase tables (in defaults)
dist = 0:1:180;
depth = 0:10:350;

eq_dist = eq_dist + 0.0001;
eq_depth = eq_depth + 0.0001;
if eq_depth>max(depth)
    eq_depth = max(depth)-0.0001;
end

% find event parameters in vectors
[~, ind_dist] = sort(abs(dist-eq_dist));
[~, ind_depth] = sort(abs(depth-eq_depth));

dist_vec(1) = dist(ind_dist(1));
dist_vec(2) = dist(ind_dist(2));

depth_vec(1) = depth(ind_depth(1));
depth_vec(2) = depth(ind_depth(2));

% find travel times in phase matrix
q1_vec(1)=phase(ind_dist(1),ind_depth(1));
q1_vec(2)=phase(ind_dist(1),ind_depth(2));

q2_vec(1)=phase(ind_dist(2),ind_depth(1));
q2_vec(2)=phase(ind_dist(2),ind_depth(2));

q3_vec(1)=interp1(depth_vec,q1_vec,eq_depth);
q3_vec(2)=interp1(depth_vec,q2_vec,eq_depth);

q_final = interp1(dist_vec,q3_vec,eq_dist);

end