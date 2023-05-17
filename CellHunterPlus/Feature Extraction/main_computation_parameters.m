clear all; close all; clc


%% tracking results loading
load(['F:\miic\roi\id20161230_pos1_roi6.mat'],'roi2save')
load(['F:\miic\roi_tu\id20161230_pos1_roi6.mat'],'roi2save_tu')
load(['F:\miic\tumor_traj rescal\tu_traj_roi_20161230_pos1_6.mat'],'track_tu')
load(['F:\miic\immune_traj rescal\sp_traj_roi_20161230_pos1_6.mat'],'track_sp')
roi = roi2save;
roi_tu = roi2save_tu; %roi cropped
% track_tu stands for MCC trajectory within the ROI under study
% track_sp stands for immune cell trajectories within the ROI under
% study
%% INTERACTION DESCRIPTORS computation
r_tu = 14; % mean MCC radius
r_sp = 4; % mean immune cell radius
r2 = 2*(r_tu+r_sp); %interaction radius r2
[count_cell_at_frame_r2, ~, sum_cum_int_r2, max_cum_int_r2,  ~, min_dist_r2, ~, ~, IDInt] = my_interaction_par_computation(track_sp,track_tu,r2);
[mean_v_r2] = my_relative_speed_sp(track_tu,track_sp,r2); 
r1 = 20; %interaction radius r1
[count_cell_at_frame_r1] = my_interaction_par_computation(track_sp,track_tu,r1);
count_cell_at_frame_r1(intersect(find(IDInt==1),find(isnan(count_cell_at_frame_r1))))=0;
[mean_v_r1] = my_relative_speed_sp(track_tu,track_sp,r1);
%% SHAPE DESCRIPTORS
Nframes = length(track_tu.t); 
m = 1;
[Area, Perimeter, Circularity, Eccentricity, InstChangeShape]= my_shape_par_computation(roi_tu,Nframes,m);
%% MOTILITY DESCRIPTORS
trajectory=[track_tu.x' track_tu.y' track_tu.t'];
[InstSpeed, StraightIndex, TotDist, NetDist]= my_motility_par_computation(trajectory(:,1:2),trajectory(:,3),Nframes);
