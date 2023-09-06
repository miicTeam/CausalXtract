% CellHunter algorithm

function [track_tu, track_sp] = cell_Hunter_Munkres(vidFrames,...
    path_temp)
% load parameters
load(path_temp, 'polarity', 'flag_imm',...
    'r_sp','Rmax_sp','DP_sp','L_sp', 'r_std',...
    'r_tu','Rmax_tu','DP_tu','dist_tu');

% Detection of cancer cells and immune cells
[data_sp,data_tu] = cell_location_ph_hunter(vidFrames,...
    r_tu, r_sp, polarity, dist_tu, flag_imm);
% Tracking algorithm based on Munkres optimization for the MCC (Main Cancer
% Cell)
[p_tu,q_tu]= link_data_Munkres(data_tu,Rmax_tu,DP_tu);
% Tracks' filtering for the MCC
[track_tu,~]= track_refining(p_tu,q_tu, 0, 0,'tu');

% Tracking algorithm and filtering for immune cells
if flag_imm
    [p_sp,q_sp]= link_data_Munkres(data_sp,Rmax_sp,DP_sp);
    [p_sp3,q_sp3]= track_refining(p_sp,q_sp,r_std,L_sp,'sp');
    track_sp=[p_sp3,q_sp3];
else
    track_sp=[];
end

end