%% Load aux functions
addpath("lib/Auxiliar/")
fem_file = "input/u_10000.h5";

%% Read info from h5 file
[t_dim, x_dim, y_dim, z_dim, u_mat_x, u_mat_z] = load_mid_u(fem_file);

%% Animation

fig = surf(x_dim, z_dim, squeeze(u_mat_z(1, :, :))');
view(30, 30)
xlabel('Distancia lateral (m)')
ylabel('Distancia axial (m)')
zlabel('Desplazamiento (m)')
zlim([min(u_mat_z, [], 'all'), max(u_mat_z, [], 'all')]*1.1)

for t_i = 1:size(u_mat_z, 1)
    fig.ZData = squeeze(u_mat_z(t_i, :, :))';
    pause(0.05)
end

%% Space-Time representation

z_plot = find(z_dim >= 10e-3, 1, 'first');

fig = figure(2);
imagesc(t_dim*1e3, x_dim*1e3, squeeze(u_mat_z(:, :, z_plot))')
colorbar
xlabel('Tiempo (ms)')
ylabel('Distancia lateral (mm)')
