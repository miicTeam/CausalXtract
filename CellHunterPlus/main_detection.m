% CellHunter+: "Segmentation and Tracking" module.
% Code for obtaining, for *each* roi,
% cancer cells' trajectories and mmune cells' trajectories;

clc, clear; close all;

% create folders for saving the following results:
% cancer cells' trajectories, immune cells' trajectories, features' data
% (computed in the "Feature Extraction" module).
path_results="../Results";
[path_tumor_traj, path_imm_traj, ~]=create_folder_result(...
    path_results);

% get rois' names
path_roi="../MATLAB_DATA/ROI MAT/";
all_names=get_files_names(path_roi);

% parameters loaded from precompiled .csv file
data = table2struct(readtable(...
    "./Config/parameters_CellHunterPlus.csv", "TextType", "string"));
check_data(data)
% parameters saved in a .mat file to easily pass them to the other module
path_name_temp=fullfile(path_results, 'temp.mat');
save(path_name_temp, '-struct', 'data');

% analyse one roi at the time
for i=1:size(all_names,1)
    n_roi=all_names{i};
    % load roi
    load(fullfile(path_roi, strcat(n_roi, '.mat')), 'roi2save'); 
    roi=roi2save;
    % segmentation and tracking
    [track_tu, track_sp] = cell_Hunter_Munkres(roi, path_name_temp);
    % save tracks
    save(fullfile(path_tumor_traj, strcat('track_tu_', n_roi, '.mat')), 'track_tu');
    if data.flag_imm
        save(fullfile(path_imm_traj, strcat('track_sp_', n_roi, '.mat')), 'track_sp');
    end

    clear n_roi roi2save roi track_tu track_sp
end
