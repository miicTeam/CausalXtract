function [p_sp2,q_sp2,data_sp,data_tu] = my_cell_Hunter_Munkres(vidFrames)
r_min_sp=4;
% Parameters setting
R_max_sp = round(5*r_min_sp);
DP_sp = 10;
lmin_sp = 10;
% Detection
[data_sp,data_tu] = my_cell_location_ph_hunter(vidFrames, r_min_sp);
% Tracking algorthm based on Munkres optimization
[p_sp,q_sp]=my_link_data_Munkres(data_sp,R_max_sp,DP_sp);
[p_sp2,q_sp2]=my_track_refining(p_sp,q_sp,2,lmin_sp,'sp');

end