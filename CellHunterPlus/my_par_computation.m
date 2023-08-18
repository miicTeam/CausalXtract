function data_matrix=my_par_computation(roi_tu, track_tu, track_sp, ...
    r2, r1, m)

%% INTERACTION DESCRIPTORS
[count_cell_at_frame_r2, ~, ~, ~,  ~, min_dist_r2, ~, ~, IDInt] = my_interaction_par_computation(track_sp,track_tu,r2);
[mean_v_r2] = my_relative_speed_sp(track_tu,track_sp,r2); 
[count_cell_at_frame_r1] = my_interaction_par_computation(track_sp,track_tu,r1);
count_cell_at_frame_r1(intersect(find(IDInt==1),find(isnan(count_cell_at_frame_r1))))=0;
[mean_v_r1] = my_relative_speed_sp(track_tu,track_sp,r1);
%% SHAPE DESCRIPTORS
Nframes = length(track_tu.t); 
[Area, Perimeter, Circularity, Eccentricity, InstChangeShape]= my_shape_par_computation(roi_tu,Nframes,m);
%% MOTILITY DESCRIPTORS
trajectory=[track_tu.x' track_tu.y' track_tu.t'];
[InstSpeed, StraightIndex, ~, NetDist]= my_motility_par_computation(trajectory(:,1:2),trajectory(:,3),Nframes);
%% Obtain matrix with features computed
data_matrix=my_feature_matrix(trajectory, ...
    Area, InstSpeed, InstChangeShape, Perimeter, Eccentricity, Circularity, NetDist,...
    StraightIndex, count_cell_at_frame_r2, min_dist_r2, mean_v_r2, mean_v_r1, count_cell_at_frame_r1);

end