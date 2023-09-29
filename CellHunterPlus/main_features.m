% CellHunter+: Feature extraction.
% Code for obtaining, for *each* MCC (Main Cancer Cell):
% interaction (with immune cells) descriptors;
% shape descriptors;
% motility descriptors
% in the form a .csv file that will be the input of tMIIC.

clc, clear; close all;

% get roi's names
path_roi="../MATLAB_DATA/ROI MAT/";
path_roi_tu='../MATLAB_DATA/ROI_TU MAT';
path_tum_traj="../results/NEW TUMOR TRAJECTORIES/";
path_imm_traj="../results/NEW IMMUNE TRAJECTORIES/";
all_names=get_files_names(path_roi);

% parameters
load('../results/temp.mat','r_tu','r_sp', 'flag_imm');
r2 = 2*(r_tu+r_sp); % interaction radius r2
r1 = r_tu+r_sp+2; % interaction radius r1

% table where the features will be saved
features_data=table();
variable_names=get_variable_names(flag_imm);

% analyse one roi at the time
for i=1:size(all_names,1)
    n_roi=all_names{i};
    % choose m based on the video
    % The video "20170517" has a different resolution than the other two, "20161230" and "20170105", 
    % m=2 compensates for that.
    if contains(n_roi, "20170517")
        m=2;
    else
        m=1;
    end
    % cropped roi and tracking results loading
    load(fullfile(path_roi_tu, strcat(n_roi, '.mat')), 'roi2save_tu')
    load(fullfile(path_tum_traj, strcat('track_tu_', n_roi, '.mat')),'track_tu')
    if flag_imm
        load(fullfile(path_imm_traj, strcat('track_sp_', n_roi, '.mat')),'track_sp')
    end
    roi_tu = roi2save_tu; % roi cropped
    % track_tu stands for MCC trajectory within the ROI under study
    % track_sp stands for immune cell trajectories within the ROI under
    % study;
    if flag_imm
        data_matrix=par_computation(roi_tu, m, track_tu, track_sp, r2, r1);
    else
        data_matrix=par_computation(roi_tu, m, track_tu);
    end
    data_matrix=[variable_names; data_matrix];
    % save result for the single roi
    writecell(data_matrix, fullfile("../Results/FEATURES/", strcat("Table_", n_roi, '.xlsx')));
    % create .csv file with all the rois' data
    features_data((1:(size(data_matrix,1)-1))+size(features_data,1),:)=data_matrix(2:end,:);
    clear n_roi roi2save_tu track_tu track_sp roi_tu m data_matrix;
end

% postprocessing
features_data.Properties.VariableNames=string(variable_names);
% add the following variables:"CAF_presence", "treatment", "apoptosis",
% "division"
file_data="./Config/state_data.csv";
if isfile(file_data)
    original_data=readtable(file_data);
    state_data=original_data(:,3:end);
    if ~isequal(all_names, unique(original_data.IDExp))
        disp("Write the state conditions in the state_data.csv file " + ...
            "following the same order used by MATLAB to upload the rois.")
        fprintf(2, '%s \n ', all_names{:})
    end
    assert(isequal(original_data(:,2), features_data(:,1)), 'There is a mistake in the order in which the ID_frame were written in the state_data.csv file.');
    input_tMIIC=[features_data(:,1), state_data, features_data(:,2:end)];
else
    input_tMIIC=features_data;
end
% save input of tMIIC
writetable(input_tMIIC, "../Demo/CausalXtract_Publication_Dataset.csv")