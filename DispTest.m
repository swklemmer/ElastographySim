addpath('lib/Auxiliar/')
addpath('lib/poly2D/')
addpath('lib/Field_II/')
rng(6942069)
field_init(0)
graf = 1;
tic()
%% Simulation parameters
% Siemens VF10-5 
% = [fs, c_c, f0, n_elem, elem_w, elv_focus, tx_focus, att, n_active] 
sim_data = [50e6, 1540, 6.67e6, 128, 0.3e-3, 20e-3, 10e-3, 70e-6, 64];

% Configurate simulation
simulation_config(sim_data);

% Configurate transducer 
trans_tx = create_transducer(sim_data);
trans_rx = create_transducer(sim_data);
sprintf("%.2f : Configuration ready", toc())

%% Create dummy displacement
real_x = (-1:0.1:10) * 1e-3;
real_y = (-1:1:1) * 1e-3;
real_z = (-1:0.1:20) * 1e-3;
[y_grid, x_grid, z_grid] = meshgrid(real_y, real_x, real_z);

u_x = zeros(length(real_x), 3, length(real_z));
u_y = zeros(length(real_x), 3, length(real_z));
disp_mag = 0.1e-3;
%u_z = disp_mag * ones(length(disp_x), 3, length(disp_z));
u_z = disp_mag * exp(-((x_grid - 3e-3).^2 + (z_grid - 3e-3).^2) / ...
      (2 * (0.5e-3)^2));

% Normalize displacement
u_z_norm = squeeze(u_z(:, 2, :)) ./ max(u_z(:, 2, :), [], 'all');

%% Iterate over parameters

%density_list = [0.1, 0.4, 0.8, 1.5, 2.5, 4, 6, 9, 12];
density_list = [6];
errors = zeros(length(density_list), 2);

for i = 1:length(density_list)
    %% Place scatterers
    % Scatterer density calculation
    cell_vol = resolution_cell(sim_data, 0);
    density = density_list(i) / (cell_vol * 1e9); % [scatters/mm^3]
    
    % Dispose scatterers at random
    phan_dim = [6, 6] * 1e-3;
    scat_pos = create_scatterers_3d([0, phan_dim(1)], [0, 0.3e-3], ...
                                    [0, phan_dim(2)], density);
    %scat_pos = [3, 0, 3; 4, 0, 4] * 1e-3;
    fprintf("\n%.2f : Scatterers ready\n", toc())

    %% Apply displacement to scatterers
    new_scat_pos = apply_u(scat_pos, real_x, real_y, real_z, u_x, u_y, u_z);
    fprintf("\n%.2f : Displacement ready\n", toc())
    
    %% Generate pre-deformation sonogram
    img_x = (0:0.05:8) * 1e-3;
    [img_z, bmode_img] = bmode_image(trans_tx, trans_rx, sim_data, ...
                                             scat_pos, img_x, 1e3);
    img_z = img_z - 0.62e-3; % Why?
    
    sonograms = zeros(2, length(img_x), length(img_z));
    sonograms(1, :, :) = bmode_img;
    fprintf("\n%.2f : Pre B-Mode ready\n", toc())

    %% Calculate post-deformation sonogram
    [~, bmode_img] = bmode_image(trans_tx, trans_rx, sim_data, ...
                                     new_scat_pos, img_x, 1e3);
    sonograms(2, :, :) = bmode_img;
    fprintf("\n%.2f : Post B-Mode ready\n", toc())

    %% Estimate displacement
    win_len = 1.5e-3;  % Window size [m]
    win_hop = 0.15e-3; % Hop between windows [m]
    fine_prec = 0.1;   % Fine precission [%/%]
    [est_x, est_z, u_x_est, u_z_est] = estimate_u(img_x, img_z, sonograms,...
                                        win_len, win_hop, fine_prec, 0);
    
    % Apply filter to estimation
    u_est_filt = conv2(...
                       medfilt2(squeeze(u_z_est(1, :, :)), [7, 7]),...
                       ones(3, 3)/9, 'same');

    % Normalize estimation
    u_est_filt = u_est_filt ./ max(u_est_filt, [], 'all');
    fprintf("\n%.2f : Disp. Estimation ready\n", toc())

    %% Calculate error after filtering
    [~, mean_err] = estimation_error(est_x, est_z, squeeze(u_z_est), ...
                                     real_x, real_z, u_z_norm, phan_dim);
    errors(i, 1) = mean_err;

    % Calculate error after filtering
    [u_z_err, mean_err] = estimation_error(est_x, est_z, u_est_filt, ...
                                     real_x, real_z, u_z_norm, phan_dim);
    errors(i, 2) = mean_err;

    %% Show results
    if ~graf
        fig = figure(1);
        fig.Position = [0, 500, 500, 300];
        show_sonogram(img_x, img_z, sonograms, 1e-23, phan_dim,...
                      1, 'First Sonogram');
    
        fig = figure(2);
        fig.Position = [300, 500, 500, 300];
        show_sonogram(img_x, img_z, sonograms, 1e-23, phan_dim,...
                      2, 'Second Sonogram');
    
        figure(3)
        fig.Position = [900, 500, 500, 300];
        show_scat_disp(scat_pos, new_scat_pos, phan_dim)
    end
    
    if graf
        fig = figure(4);
        fig.Position = [0, 100, 1500, 500];
        [~, ~, ~] = show_disp_error(real_x, real_z, u_z_norm, ...
                                    est_x, est_z, u_est_filt,...
                                    u_z_err, errors(i, 2), phan_dim);
        %saveas(fig, sprintf("out_fig/est_d%.1f_win_%.2f.fig", density, win_len*1e3))
    end
end

save('out_fig/density_experiment.mat', 'errors', 'density_list')


%% Terminate program and free memory 
xdc_free(trans_tx)
xdc_free(trans_rx)
field_end()
