% Correct the cancer cell's trajectory in order to decrease flickering.

function final_traj=track_correction(roi2save, path_name_temp, track_tu, indexes)

% load parameters
load(path_name_temp,...
    'DP_tu', 'Rmax_tu', 'r_tu', 'polarity', 'dist_tu');
% Detection and tracking
[data_tu]=mitosis_cell_location_ph_hunter(roi2save, indexes, r_tu, polarity, dist_tu);
[p,~]=link_data_Munkres(data_tu,Rmax_tu,DP_tu);
% Track correction
field_names=fieldnames(track_tu);
final_traj=struct();
for j=1:size(field_names,1)
    field=field_names{j};
    final_traj.(field)(1:indexes(1)-1)=track_tu.(field)(1:indexes(1)-1);
    final_traj.(field)(indexes(1):indexes(end))=p.(field)(indexes(1):indexes(end));
    clear field;
end

end
