%% Code for obtaining all the names of the files (rois) in the folder.
function [all_names, dates, posits, rois]=my_files_names(path_roi)
all_files=dir(fullfile(path_roi, "*.mat"));
all_names = {all_files.name}';
dates=cellfun(@(x) extractBetween(x,'id','_'), all_names, 'UniformOutput', true);
posits=cellfun(@(x) extractBetween(x,'pos','_'), all_names, 'UniformOutput', true);
rois=cellfun(@(x) extractBetween(x,'roi','.mat'), all_names, 'UniformOutput', true);
end