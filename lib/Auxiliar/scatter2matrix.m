function [x_dim, z_dim, u_matrix] = scatter2matrix(x_scatter, z_scatter, ...
                                                    u_scatter)
%SCATTER2MATRIX Transforms displacement in scatter vector form to 3D matrix
%representation.

x_dim = unique(x_scatter);
z_dim = unique(z_scatter);

dx = x_dim(2) - x_dim(1);
dz = z_dim(2) - z_dim(1);

u_matrix = zeros(size(u_scatter, 1), length(x_dim), length(z_dim));

for t = 1:size(u_scatter, 1)
    for i = 1:length(u_scatter)
        x_ind = round(1 + x_scatter(i) / dx);
        z_ind = round(1 + z_scatter(i) / dz);
        u_matrix(t, x_ind, z_ind) = u_scatter(t, i);
    end
end
end

