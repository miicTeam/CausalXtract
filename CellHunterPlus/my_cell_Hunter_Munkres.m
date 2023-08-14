%% CellHunter algorithm
function [p_sp,q_sp,p_tu,q_tu,data_sp,data_tu] = my_cell_Hunter_Munkres(vidFrames,...
    DP_sp,R_max_sp,R_max_tu,DP_tu,...
    r_min_tu_int,r_min_sp,polarity, dist_tu)
% Detection
[data_sp,data_tu] = my_cell_location_ph_hunter(vidFrames,...
    r_min_tu_int, r_min_sp, polarity, dist_tu);
% Tracking algorthm based on Munkres optimization
[p_sp,q_sp]=my_link_data_Munkres(data_sp,R_max_sp,DP_sp);
[p_tu,q_tu]=my_link_data_Munkres(data_tu,R_max_tu,DP_tu);
end