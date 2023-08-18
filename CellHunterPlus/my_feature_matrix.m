function data_matrix=my_feature_matrix(...
    trajectory, ...
    Area, InstSpeed, InstChangeShape, Perimeter, Eccentricity, Circularity, NetDist,...
    StraightIndex, count_cell_at_frame_r2, min_dist_r2, mean_v_r2, mean_v_r1, count_cell_at_frame_r1)

% change NaN to 0 when necessary
count_cell_at_frame_r2(isnan(count_cell_at_frame_r2))=0;
count_cell_at_frame_r1(isnan(count_cell_at_frame_r1))=0;
if isnan(NetDist(1))
    NetDist(1)=0;
end
if isnan(StraightIndex(1))
    StraightIndex(1)=0;
end

N_frames=trajectory(:,3);
data_matrix=[N_frames, Area', InstSpeed', InstChangeShape, Perimeter', Eccentricity', Circularity', NetDist',...
    StraightIndex', count_cell_at_frame_r2', min_dist_r2', mean_v_r2', mean_v_r1', count_cell_at_frame_r1'];
% convert NaN to NA for tMIIC
data_matrix=num2cell(data_matrix);
data_matrix(cellfun(@isnan,data_matrix)) = {"NA"};

end
