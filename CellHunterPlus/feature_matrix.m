% Return the feature matrix for each roi.

function data_matrix=feature_matrix(...
    trajectory, ...
    Area, InstSpeed, InstChangeShape, Perimeter, Eccentricity, Circularity, NetDist,...
    StraightIndex, varargin)

p = inputParser();
addOptional(p, 'count_cell_at_frame_r2', NaN, @isnumeric) 
addOptional(p, 'min_dist_r2', NaN, @isnumeric) 
addOptional(p, 'mean_v_r2', NaN, @isnumeric) 
addOptional(p, 'mean_v_r1', NaN, @isnumeric) 
addOptional(p, 'count_cell_at_frame_r1', NaN, @isnumeric) 
parse(p, varargin{:})
ip = p.Results;

% change NaN to 0 when necessary
if isnan(NetDist(1))
    NetDist(1)=0;
end
if isnan(StraightIndex(1))
    StraightIndex(1)=0;
end

N_frames=trajectory(:,3);

if size(varargin,2)==5
    % change NaN to 0 when necessary
    ip.count_cell_at_frame_r2(isnan(ip.count_cell_at_frame_r2))=0;
    ip.count_cell_at_frame_r1(isnan(ip.count_cell_at_frame_r1))=0;
    data_matrix=[N_frames, Area', InstSpeed', InstChangeShape, Perimeter', Eccentricity', Circularity', NetDist',...
        StraightIndex', ip.count_cell_at_frame_r2', ip.min_dist_r2', ip.mean_v_r2', ip.mean_v_r1', ip.count_cell_at_frame_r1'];
else
    data_matrix=[N_frames, Area', InstSpeed', InstChangeShape, Perimeter', Eccentricity', Circularity', NetDist',...
        StraightIndex'];
end

% convert NaN to NA for tMIIC
data_matrix=num2cell(data_matrix);
data_matrix(cellfun(@isnan,data_matrix)) = {"NA"};

end
