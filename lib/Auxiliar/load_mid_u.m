function [t_dim, x_dim, y_dim, z_dim, u_mat_x, u_mat_z] = load_mid_u(fem_file)
%LOAD_U Loads displacements in the middle slice from h5 file.

% Read node positions from h5 file
node_positions = h5read(fem_file, "/Mesh/mesh/geometry");

% Read adquisition times from h5 file
fem_info = h5info(fem_file, "/Function/Displacement");

adq_times = struct2cell(fem_info.Datasets);
adq_times = adq_times(1, :);

% Convert times to floats and sort them
[t_dim, sort_ind] = time2float(adq_times);

% Read node z-displacements from h5 file
u_x = zeros(length(adq_times), size(node_positions, 2));
u_z = zeros(length(adq_times), size(node_positions, 2));

for i = 1:length(adq_times)
    vec_u = h5read(fem_file, strcat("/Function/Displacement/", ...
                    adq_times(i)));
    u_x(i, :) = vec_u(1, :);
    u_z(i, :) = vec_u(3, :);
end

% Find nodes in middle slice
middle_slice = node_positions(2, :) == 0;
mid_pts = node_positions([1 3], middle_slice);

% Convert scatter data into matrix form
[~, ~, u_mat_x] = scatter2matrix(mid_pts(1, :), mid_pts(2, :), ...
    u_x(sort_ind, middle_slice));

[x_dim, z_dim, u_mat_z] = scatter2matrix(mid_pts(1, :), mid_pts(2, :), ...
    u_z(sort_ind, middle_slice));

% Find Y-dimension
y_dim = unique(node_positions(2, :));
end

