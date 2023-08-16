%% CellHunter+: "Segmentation and Tracking" module.
% Code for obtaining, for *each* roi:
% immune cells' trajectories;
% cancer cells' trajectories.

clc, clear; close all;

% create folders for saving the following results:
% cancer cells' trajectories, immune cells' trajectories, features' data
% (computed in the "Feature Extraction" module).
path_results="../results";
[path_tumor_traj, path_imm_traj, ~]=my_folder_result(...
    path_results);

% get rois' names
path_roi="../MATLAB_DATA/ROI MAT/";
[all_names, dates, posits, rois]=my_files_names(path_roi);

% parameters loaded from precompiled .csv file
[polarity, r_sp, Rmax_sp, DP_sp, L_sp, r_std,...
    r_tu, Rmax_tu, DP_tu, dist_tu]= readvars(...
    "../MATLAB_DATA/parameters_CellHunterPlus.csv", "TextType", "string");
%  parameters saved in .mat file to easily pass the parameters to the other
%  module
save(fullfile(path_results, "temp.mat"),'polarity',...
    'r_sp','Rmax_sp','DP_sp','L_sp', 'r_std',...
    'r_tu','Rmax_tu','DP_tu','dist_tu');

% analyse one roi at the time
for i=1:size(all_names,1)
    dt=dates{i};
    n_pos=posits{i};
    n_roi=rois{i};
    % load roi
    load(fullfile(path_roi, sprintf('id%s_pos%s_roi%s.mat', dt, n_pos, n_roi)), 'roi2save'); 
    roi=roi2save;
    % segmentation and tracking
    [p_sp2,q_sp2,p_tu2,q_tu2,~,~] = my_cell_Hunter_Munkres(roi,...
        DP_sp,Rmax_sp,Rmax_tu,DP_tu,r_tu,r_sp,polarity, dist_tu);
    % tracks' filtering
    [p_sp3,q_sp3]=my_track_refining(p_sp2,q_sp2,r_std,L_sp,'sp');
    track_sp=[p_sp3,q_sp3];
    [track_tu,~]=my_track_refining(p_tu2,q_tu2, 0, 0,'tu');
    % save tracks
    save(fullfile(path_imm_traj, sprintf("new_sp_traj_roi_%s_pos%s_%s.mat", dt, n_pos, n_roi)), 'track_sp');
    save(fullfile(path_tumor_traj, sprintf("new_tu_traj_roi_%s_pos%s_%s.mat", dt, n_pos, n_roi)), 'track_tu');
    
    clear dt n_pos n_roi roi2save roi...
        p_sp2 q_sp2 p_tu2 q_tu2...
        p_sp3 q_sp3 track_sp track_tu
end
