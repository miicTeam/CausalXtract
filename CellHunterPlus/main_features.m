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
load('../results/temp.mat','r_tu','r_sp', 'flag_imm');
r2 = 2*(r_tu+r_sp); % interaction radius r2
r1 = r_tu+r_sp+2; % interaction radius r1

% table where the features will be saved
features_data=table();
variable_names=my_variable_names(flag_imm);

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
    if flag_imm
        load(sprintf('../results/NEW IMMUNE TRAJECTORIES/new_sp_traj_roi_%s_pos%s_%s.mat', dt, n_pos, n_roi),'track_sp')
    end
    roi_tu = roi2save_tu; % roi cropped
    % track_tu stands for MCC trajectory within the ROI under study
    % track_sp stands for immune cell trajectories within the ROI under
    % study;
    if flag_imm
        data_matrix=my_par_computation(roi_tu, m, track_tu, track_sp, r2, r1);
    else
        data_matrix=my_par_computation(roi_tu, m, track_tu);
    end
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
file_data="./Config/state_data.csv";
if isfile(file_data)
    original_data=readtable(file_data);
    state_data=original_data(:,3:end);
    if ~isequal(cellfun(@(x) x(1:end-4), all_names, 'UniformOutput', false), unique(original_data.IDExp))
        disp("Write the state conditions in the state_data.csv file following the same order used by MATLAB to upload the rois.")
        fprintf(2, '%s \n ', all_names{:})
    end
    assert(isequal(original_data(:,2), features_data(:,1)), 'There is a mistake in the order in which the ID_frame were written in the state_data.csv file.');
    input_tMIIC=[features_data(:,1), state_data, features_data(:,2:end)];
else
    input_tMIIC=features_data;
end
% save input of tMIIC
writetable(input_tMIIC, "../Demo/CausalXtract_Publication_Dataset.csv")