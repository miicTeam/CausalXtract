% Obtain all the names of the files (rois) in the folder.

function roi_names=get_files_names(path_roi)
all_files=dir(fullfile(path_roi, "*.mat"));
file_names = {all_files.name}';
roi_names=cellfun(@(x) extractBefore(x, '.mat'), file_names, 'UniformOutput', false);
end