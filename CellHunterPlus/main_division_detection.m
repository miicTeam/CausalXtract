% Code for refining cell's detection and tracking in case of mitosis.
% The objective is to:
% identify both daughter cells originating from a mother cell;
% define the position of a fictitious cell, whose coordinates (x, y) are
% the mean of both daughter cells' coordinates.
% This allows to decrease the "flickering" that characterizes the MCC's
% trajectory when it divides, if the two daughter cells are correctly
% identified.
clc; clear; close all
% extract names of files
path_roi="../MATLAB_DATA/ROI MAT/";
[all_names, dates, posits, rois]=my_files_names(path_roi);
path_tum_traj="../results/NEW TUMOR TRAJECTORIES/";
% load division data
T=readtable("../MATLAB_DATA/global_division.csv");
glob_division=T.global_division;
frame_division=T.frame_division;
% load parameters
data=load('../results/temp.mat','DP_tu','Rmax_tu', 'dist_tu');
DP_tu=data.DP_tu;
Rmax_tu=data.Rmax_tu;
dist_tu=data.dist_tu;
% apply correction
for i=1:size(all_names,1)
    if glob_division(i)==1
        dt=dates{i};
        n_pos=posits{i};
        n_roi=rois{i};
        % load roi and cancer cell trajectory
        load(fullfile(path_roi, sprintf('id%s_pos%s_roi%s.mat', dt, n_pos, n_roi)), 'roi2save');
        load(fullfile(path_tum_traj, sprintf('new_tu_traj_roi_%s_pos%s_%s.mat', dt, n_pos, n_roi)),'track_tu')
        % correct after mitosis
        indexes=frame_division(i):size(roi2save,3);
        [data_tu]=my_mitosis_cell_location_ph_hunter(roi2save, indexes, dist_tu);
        [p,~]=my_link_data_Munkres(data_tu,Rmax_tu,DP_tu);
        field_names=fieldnames(track_tu);
        final_traj=struct();
        for j=1:size(field_names,1)
            field=field_names{j};
            final_traj.(field)(1:indexes(1)-1)=track_tu.(field)(1:indexes(1)-1);
            final_traj.(field)(indexes(1):indexes(end))=p.(field)(indexes(1):indexes(end));
            clear field;
        end
        clear track_tu;
        track_tu=final_traj;
        save(fullfile(path_tum_traj, sprintf("new_tu_traj_roi_%s_pos%s_%s.mat", dt, n_pos, n_roi)), 'track_tu');
    end
    clear dt n_roi n_pos...
        roi2save track_tu...
        indexes data_tu p field_names final_traj;
end
