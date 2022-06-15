function scat_pos = create_scatterers_3d(x_dim, y_dim, z_dim, density)
%CREATE_SCATTERERS Dispose scatterers at random in the region established
%by given dimension. Asumes plane symmetry in the x and y dimension.

scat_pos = rand(floor(density * 1e9 * ...
                x_dim(end) * 2 * y_dim(end) * z_dim(end)), 3);

scat_pos = scat_pos * diag([x_dim(end), 2 * y_dim(end), z_dim(end)]);

% Center phantom
scat_pos = scat_pos + [0, -y_dim(end), 2e-3];

end

