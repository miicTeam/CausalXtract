% Code for refining cell's detection and tracking in case of mitosis.
% The objective is to:
% identify both daughter cells originating from a mother cell;
% define the position of a fictitious cell, whose coordinates (x, y) are
% the mean of both daughter cells' coordinates.
% This allows to decrease the "flickering" that characterizes the MCC's
% trajectory when it divides, if the two daughter cells are correctly
% identified.

clc; clear; close all

% load division data
name_file="./Config/global_division.csv";
if isfile(name_file)
    T=readtable(name_file);
    IDExp=T.IDExp;
    glob_division=T.global_division;
    frame_division=T.frame_division;

    % extract names of files
    path_roi="../MATLAB_DATA/ROI MAT/";
    all_names=my_files_names(path_roi);
    path_tum_traj="../Results/NEW TUMOR TRAJECTORIES/";
    path_name_temp="../Results/temp.mat";
    
    % apply correction
    for i=1:size(all_names,1)
        idx=contains(all_names,IDExp{i});
        if glob_division(idx)==1
            n_roi=all_names{idx};
            % load roi and cancer cell trajectory
            load(fullfile(path_roi, strcat(n_roi, '.mat')), 'roi2save');
            load(fullfile(path_tum_traj, strcat('track_tu_', n_roi, '.mat')),'track_tu')
            roi=roi2save;
            % correct after mitosis
            indexes=frame_division(idx):size(roi,3);
            [final_traj]=my_track_correction(roi, path_name_temp, track_tu, indexes);
            clear track_tu;
            track_tu=final_traj;
            save(fullfile(path_tum_traj, strcat('track_tu_', n_roi, '.mat')), 'track_tu');
        end
        clear n_roi roi2save roi track_tu indexes final_traj;
    end
end