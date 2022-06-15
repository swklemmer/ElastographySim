function new_scat_pos = apply_u(scat_pos, x_dim, y_dim, z_dim, ...
                                u_x, u_y, u_z)
%APPLY_U modify scatterers position using the displacement output from the
%FEM simulation.

% Substract safety gap
scat_pos = scat_pos - [0, 0, 2e-3];
new_scat_pos = scat_pos;

% Interpolate values into new scatterer array
new_scat_pos(:, 1) = scat_pos(:, 1) + ...
        interp3(y_dim, x_dim, z_dim, u_x, abs(scat_pos(:, 2)), ...
        abs(scat_pos(:, 1)), scat_pos(:, 3), "linear");

new_scat_pos(:, 2) = scat_pos(:, 2) + ...
        interp3(y_dim, x_dim, z_dim, u_y, abs(scat_pos(:, 2)), ...
        abs(scat_pos(:, 1)), scat_pos(:, 3), "linear");

new_scat_pos(:, 3) = scat_pos(:, 3) + ...
        interp3(y_dim, x_dim, z_dim, u_z, abs(scat_pos(:, 2)), ...
        abs(scat_pos(:, 1)), scat_pos(:, 3), "linear");

% Add safety gap
new_scat_pos = new_scat_pos + [0, 0, 2e-3];

end

