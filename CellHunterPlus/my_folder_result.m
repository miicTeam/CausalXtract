% Create folders for saving the results of CellHunter+.
% One folder for tumor trajectories, one for immune trajectories, one for
% features' data (feature extraction module).

function [path_tumor_traj, path_imm_traj, path_tables]=my_folder_result(...
    path_results)
path_tumor_traj=fullfile(path_results, "NEW TUMOR TRAJECTORIES");
path_imm_traj=fullfile(path_results, "NEW IMMUNE TRAJECTORIES");
path_tables=fullfile(path_results, "FEATURES");
if not(isfolder(path_results))
    mkdir(path_results);
end
if not(isfolder(path_tumor_traj))
    mkdir(path_tumor_traj);
end
if not(isfolder(path_imm_traj))
    mkdir(path_imm_traj);
end
if not(isfolder(path_tables))
    mkdir(path_tables);
end
end