% file of example to apply tracking to a crop
clear; close all; clc;
info=imfinfo('E:\CASUAL_X_TRACK\MIIC_OLD_DATA\id20161230_pos1_roi6.tif');
for k = 1 : length(info)
    video(:,:,k)=im2double(imread(info(k).Filename,'Index',k));
end

[p_sp2,q_sp2,data_sp,data_tu] = my_cell_Hunter_Munkres(video);