function new_scat_pos = apply_mid_u(scat_pos, x_dim, z_dim, u_z)
%APPLY_U modify scatterers position using the displacement output from the
%FEM simulation.

% Create mesh grid
[x_grid, z_grid] = meshgrid(x_dim, z_dim);

% Interpolate values into new scatterer array
new_scat_pos = scat_pos;

new_scat_pos(:, 3) = scat_pos(:, 3) + ...
        interp2(x_grid, z_grid, u_z', abs(scat_pos(:, 1)), scat_pos(:, 3), ...
        "linear");

end

