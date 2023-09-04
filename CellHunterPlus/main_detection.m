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
data = table2struct(readtable(...
    "./Config/parameters_CellHunterPlus.csv", "TextType", "string"));
check_data(data)
% parameters saved in a .mat file to easily pass them to the other module
path_name_temp=fullfile(path_results, 'temp.mat');
save(path_name_temp, '-struct', 'data');

% analyse one roi at the time
for i=1:size(all_names,1)
    dt=dates{i};
    n_pos=posits{i};
    n_roi=rois{i};
    % load roi
    load(fullfile(path_roi, sprintf('id%s_pos%s_roi%s.mat', dt, n_pos, n_roi)), 'roi2save'); 
    roi=roi2save;
    % segmentation and tracking
    [track_tu, track_sp] = my_cell_Hunter_Munkres(roi, path_name_temp);
    % save tracks
    save(fullfile(path_tumor_traj, sprintf("new_tu_traj_roi_%s_pos%s_%s.mat", dt, n_pos, n_roi)), 'track_tu');
    if data.flag_imm
        save(fullfile(path_imm_traj, sprintf("new_sp_traj_roi_%s_pos%s_%s.mat", dt, n_pos, n_roi)), 'track_sp');
    end

    clear dt n_pos n_roi roi2save roi track_sp track_tu
end
