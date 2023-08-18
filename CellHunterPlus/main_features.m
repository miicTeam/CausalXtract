%% CellHunter+: Feature extraction.
% Code for obtaining, for *each* MCC (Main Cancer Cell):
% interaction (with immune cells) descriptors;
% shape descriptors;
% motility descriptors
% in the form a .csv file that will be the input of tMIIC.

clc, clear; close all;

% get roi's names
path_roi="../MATLAB_DATA/ROI MAT/";
[all_names, dates, posits, rois]=my_files_names(path_roi);

% parameters
data=load('../results/temp.mat','r_tu','r_sp');
r_tu = data.r_tu; % mean MCC radius
r_sp = data.r_sp; % mean immune cell radius
r2 = 2*(r_tu+r_sp); % interaction radius r2
r1 = r_tu+r_sp+2; % interaction radius r1

% table where the features will be saved
features_data=table();
variable_names={"ID_frame",...
    "area", "instantaneous_cancer_velocity",...
    "instantaneous_shape_change", "perimeter", "eccentricity", "circularity",...
    "net_displacement", "directionality", "number_of_cancer_immune_interactions_r2",...
    "minimal_cancer_immune_distance_r2", "mean_immune_velocity_r2",...
    "mean_immune_velocity_r1", "number_of_cancer_immune_interactions_r1"};

% analyse one roi at the time
for i=1:size(all_names,1)
    dt=dates{i};
    n_pos=posits{i};
    n_roi=rois{i};
    % choose m based on the video
    if strcmp(dt, "20170517")
        m=2;
    elseif or(strcmp(dt, "20161230"), strcmp(dt, "20170105"))
        m=1;
    end
    % roi and tracking results loading
    load(sprintf('../MATLAB_DATA/ROI_TU MAT/id%s_pos%s_roi%s.mat', dt, n_pos, n_roi),'roi2save_tu')
    load(sprintf('../results/NEW TUMOR TRAJECTORIES/new_tu_traj_roi_%s_pos%s_%s.mat', dt, n_pos, n_roi),'track_tu')
    load(sprintf('../results/NEW IMMUNE TRAJECTORIES/new_sp_traj_roi_%s_pos%s_%s.mat', dt, n_pos, n_roi),'track_sp')
    roi_tu = roi2save_tu; % roi cropped
    % track_tu stands for MCC trajectory within the ROI under study
    % track_sp stands for immune cell trajectories within the ROI under study
    data_matrix=my_par_computation(roi_tu, track_tu, track_sp, r2, r1, m);
    data_matrix=[variable_names; data_matrix];
    % save result for the single roi
    writecell(data_matrix, sprintf("../results/tables/table_%s_pos%s_%s.xlsx", dt, n_pos, n_roi)) 
    % create .csv file with all the rois' data
    features_data((1:(size(data_matrix,1)-1))+size(features_data,1),:)=data_matrix(2:end,:);
    clear dt n_pos n_roi roi2save_tu track_tu track_sp roi_tu m data_matrix;
end

% postprocessing
features_data.Properties.VariableNames=string(variable_names);
% add the following variables:"CAF_presence", "treatment", "apoptosis",
% "division"
original_data=readtable("../MATLAB_DATA/state_data.csv");
state_data=original_data(:,2:end);
input_tMIIC=[features_data(:,1), state_data, features_data(:,2:end)];
% save input of tMIIC
writetable(input_tMIIC, "../Demo/CausalXtract_Publication_Dataset.csv")