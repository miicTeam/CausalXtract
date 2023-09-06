% Check that the data is correctly inserted in the .csv file.

function check_data(data)
% check if the required parameter are correctly inserted in the .csv file.
name_base_vars=["polarity", "flag_imm", "r_tu", "Rmax_tu", "DP_tu", "dist_tu"];
check_presence_data(data, name_base_vars);

if ~any(strcmp(data.polarity, ["bright", "dark"]))
    error('Insert parameter: polarity, either as "bright" or "dark".')
end

if ~or(data.flag_imm==0, data.flag_imm==1)
    error('Insert parameter: flag_imm, either as 0 or 1.')
end

for j=3:length(name_base_vars)
    validateattributes(data.(name_base_vars(j)), "numeric", {"nonempty", "positive"})
end

% check if immune variables parameters are correctly inserted in the file
if data.flag_imm
    name_imm_vars=["r_sp", "Rmax_sp", "DP_sp", "L_sp", "r_std"];
    check_presence_data(data, name_imm_vars);
    for i=1:length(name_imm_vars)
        validateattributes(data.(name_imm_vars(i)), "numeric", {"nonempty", "positive"})
    end
end


function check_presence_data(data, name_vars)
struct_vars=struct();
for i=1:length(name_vars)
    struct_vars.(name_vars(i))=isfield(data, name_vars(i));
    if ~(struct_vars.(name_vars(i)))
        error("Please insert the following parameter in the .csv file: %s", name_vars(i))
    end
end

