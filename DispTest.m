addpath('lib/Auxiliar/')
addpath('lib/poly2D/')
addpath('lib/Field_II/')
field_init(0)
graf = 0;
tic()
%% Simulation parameters
% Siemens VF10-5 = [fs, c_c, f0, n_elem, elem_w, elv_focus, tx_focus, alpha] 
sim_data = [50e6, 1540, 6.67e6, 128, 0.3e-3, 10e-3, 10e-3, 70e-6];

% Configurate simulation
simulation_config(sim_data);

% Configurate transducer 
trans_tx = create_transducer(sim_data);
trans_rx = create_transducer(sim_data);
sprintf("%.2f : Configuration ready", toc())

%% Place scatterers
% Dispose scatterers at random
density = 5; % [scatters/mm^3]
scat_pos = create_scatterers_3d([0, 6e-3], [0, 0.5e-3], [0, 10e-3], density);

sprintf("%.2f : Scatterers ready", toc())

%% Generate pre-deformation sonogram
img_x = (0:199) * 0.5e-4;
n_active = 63;
[img_z, bmode_img] = bmode_image(trans_tx, trans_rx, sim_data, ...
                                         scat_pos, img_x, 1.5e3, n_active);
img_z = img_z - 1e-3;

sonograms = zeros(2, length(img_x), length(img_z));
sonograms(1, :, :) = bmode_img;
sprintf("%.2f : Pre B-Mode ready", toc())

%% Create dummy displacement
[y_grid, x_grid, z_grid] = meshgrid([-1, 0, 1] * 1e-3, img_x, img_z);

u_x = zeros(length(img_x), 3, length(img_z));
u_y = zeros(length(img_x), 3, length(img_z));
disp_mag = 0.5e-3;
%u_z = disp_mag * ones(length(img_x), 3, length(img_z));
u_z = disp_mag * exp(-((x_grid - 3e-3).^2 + (z_grid - 5e-3).^2) / (2 * (1e-3)^2));

%% Apply displacement to scatterers
new_scat_pos = apply_u(scat_pos, img_x, [-1 0 1] * 1e-3, img_z, u_x, u_y, u_z);
sprintf("%.2f : Displacement ready", toc())

%% Calculate post-deformation sonogram
[~, bmode_img] = bmode_image(trans_tx, trans_rx, sim_data, ...
                                 new_scat_pos, img_x, 1.5e3, n_active);
sonograms(2, :, :) = bmode_img;

% Terminate program and free memory 
xdc_free(trans_tx)
xdc_free(trans_rx)
field_end()
sprintf("%.2f : Post B-Mode ready", toc())

%% Estimate displacement
fine_prec = 0.1;
[est_x, est_z, u_x_est, u_z_est] = estimate_u(img_x, img_z, sonograms,...
                                              fine_prec, graf);
sprintf("%.2f : Disp. Estimation ready", toc())

%% Show results

if ~graf
    fig = figure(1);
    fig.Position = [0, 0, 500, 500];
    imagesc(img_z, img_x, squeeze(sonograms(1, :, :)));
    hold on
    scatter(scat_pos(:, 3), scat_pos(:, 1), [],'red')
    hold off
    ylabel('Lateral Distance (X) [m]')
    xlabel('Depth (Z) [m]')
    title('First Sonogram')
    colormap(jet)
    colorbar

    fig = figure(2);
    fig.Position = [500, 0, 500, 500];
    imagesc(img_z, img_x, squeeze(sonograms(2, :, :)));
    hold on
    scatter(new_scat_pos(:, 3), new_scat_pos(:, 1), [],'red')
    hold off
    ylabel('Lateral Distance (X) [m]')
    xlabel('Depth (Z) [m]')
    title('Second Sonogram')
    colormap(jet)
    colorbar

    figure(3)
    scatter(scat_pos(:, 3), scat_pos(:, 1))
    hold on
    scatter(new_scat_pos(:, 3), new_scat_pos(:, 1))
    hold off
    axis ij
    axis([img_z(1), img_z(end), img_x(1), img_x(end)])
    title('Scatterer Position')
    xlabel('Depth (m)')
    ylabel('Lateral position (m)')
    legend({'Pre', 'Post'})
end

if graf
    fig = figure(4);
    fig.Position = [0, 0, 500, 500];
    fig_mesh1 = mesh(est_z, est_x, medfilt2(squeeze(u_z_est(1, :, :)), [9, 9]));
    zlim([0, 5 * disp_mag])
    clim([0, 5 * disp_mag])
    title('Estimated Z Displacement')
    view(30, 30)
    colormap("jet")
    colorbar
    
    fig = figure(5);
    fig.Position = [500, 0, 500, 500];
    fig_mesh2 = mesh(squeeze(z_grid(:, 2, :)), squeeze(x_grid(:, 2, :)), squeeze(u_z(:, 2, :)));
    zlim([0, 5 * disp_mag])
    clim([0, 5 * disp_mag])
    title('Real Z Displacement')
    view(30, 30)
    colormap("jet")
    colorbar
    
    
    for t_i = 2:size(u_z_est, 1)
        break
        pause(0.1)
        fig_mesh1.ZData = squeeze(u_z_est(t_i, :, :))';
        fig_mesh2.ZData = squeeze(u_z(t_i, :, :))';
    end
end
