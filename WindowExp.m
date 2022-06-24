addpath('lib/Auxiliar/')
addpath('lib/poly2D/')
addpath('lib/Field_II/')
rng(6942069)
field_init(0)

%% Simulation parameters
% Siemens VF10-5 
% = [fs, c_c, f0, n_elem, elem_w, elv_focus, tx_focus, att, n_active] 
sim_data = [50e6, 1540, 6.67e6, 128, 0.3e-3, 20e-3, 10e-3, 70e-6, 64];

% Configurate simulation
simulation_config(sim_data);

% Configurate transducer 
trans_tx = create_transducer(sim_data);
trans_rx = create_transducer(sim_data);

% Scatterer density calculation
cell_vol = resolution_cell(sim_data, 0);
density = 4 / (cell_vol * 1e9); % [scatterers/mm^3]
phan_dim = [6, 6] * 1e-3;

%% Create real displacement

 % Generate displacement grid
real_x = (-1:0.1:10) * 1e-3;
real_y = (-1:1:1) * 1e-3;
real_z = (-1:0.1:20) * 1e-3;
[y_grid, x_grid, z_grid] = meshgrid(real_y, real_x, real_z);

% Write real displacement
disp_mag = 0.1e-3;

u_x = zeros(length(real_x), 3, length(real_z));
u_y = zeros(length(real_x), 3, length(real_z));
u_z = disp_mag * exp(-((x_grid - 3e-3).^2 + (z_grid - 3e-3).^2) / ...
      (2 * (0.5e-3)^2));

% Normalize displacement
u_z_norm = squeeze(u_z(:, 2, :)) ./ max(u_z(:, 2, :), [], 'all');

% Calculate Null estimation error
[~, null_err] = estimation_error(real_x, real_z, u_z_norm, ...
                                 real_x, real_z, 0 * u_z_norm, ...
                                 phan_dim);

%% Iterate over window size
win_list = (0.5:0.5:2) * 1e-3;
errors = zeros(2, length(win_list), 20);

%%
for j = 1:20

% Dispose scatterers at random
scat_pos = create_scatterers_3d([0, phan_dim(1)], [0, 0.3e-3], ...
                                [0, phan_dim(2)], density);

% Apply displacement to scatterers
new_scat_pos = apply_u(scat_pos, real_x, real_y, real_z, u_x, u_y, u_z);

% Generate pre-deformation sonogram
img_x = (0:0.1:8) * 1e-3;
[img_z, bmode_img] = bmode_image(trans_tx, trans_rx, sim_data, ...
                                         scat_pos, img_x, 1.5e3);
img_z = img_z - 0.62e-3; % WHY?

sonograms = zeros(2, length(img_x), length(img_z));
sonograms(1, :, :) = bmode_img;

% Calculate post-deformation sonogram
[~, bmode_img] = bmode_image(trans_tx, trans_rx, sim_data, ...
                                 new_scat_pos, img_x, 1.5e3);
sonograms(2, :, :) = bmode_img;

for i = 1:length(win_list)
    tic()
    win_len = win_list(i);  % Window size [m]
    win_hop = 0.15e-3; % Hop between windows [m]
    fine_prec = 0.1;   % Fine precission [%/%]

    % Estimate displacement
    [est_x, est_z, ~, u_z_est] = estimate_u(img_x, img_z, sonograms, ...
                                        win_len, win_hop, fine_prec, 0);

    % Normalize estimation
    u_z_est = u_z_est ./ max(u_z_est, [], 'all');

    % Apply filter to estimation
    u_est_filt = conv2(medfilt2(squeeze(u_z_est(1, :, :)), [7, 7]),...
                       ones(3, 3)/9, 'same');

    % Calculate error before filtering
    [~, mean_err] = estimation_error(real_x, real_z, u_z_norm, ...
                                     est_x, est_z, squeeze(u_z_est), ...
                                     phan_dim);
    errors(1, i, j) = mean_err / null_err;

    % Calculate error after filtering
    [~, mean_err] = estimation_error(real_x, real_z, u_z_norm, ...
                                     est_x, est_z, u_est_filt, ...
                                     phan_dim);
    errors(2, i, j) = mean_err / null_err;

    fprintf("\nt = %.2f: win.len. = %.2f mm: iter. = %d\n", ...
            toc(), 1e3 * win_list(i), j)
end
end

save(sprintf("out_experiments/random2_win_experiment_d%.1f.mat", density),...
    'errors', 'win_list')

%% Terminate program and free memory 
xdc_free(trans_tx)
xdc_free(trans_rx)
field_end()

%% Show results
load(sprintf("out_experiments/random_win_experiment_d%.1f.mat", density),...
    'errors', 'win_list')

fig = figure(1);
fig.Position = [500, 500, 500, 300];
boxchart(squeeze(errors(2, :, :))' / null_err);

ylim([0, 5])
set(gca,'XTickLabel', win_list * 1e3);
grid on
ylabel("Mean Absolute Error / Null Error")
xlabel("Window length [mm]")
title("Effect of window length on displacement estimation error")

%%
load(sprintf("out_experiments/random2_win_experiment_d%.1f.mat", density),...
    'errors', 'win_list')

error_vec = [reshape(errors(1, :, :), [], 1); reshape(errors(2, :, :), [], 1)];
error_win = repmat((1:length(win_list))', 40, 1);
error_filt = [zeros(length(error_vec)/2, 1); ones(length(error_vec)/2, 1)];

fig = figure(2);
fig.Position = [500, 500, 500, 300];
boxchart(error_win, error_vec, 'GroupByColor', error_filt);

ylim([0, 2])
legend({'Pre-Filter', 'Post-Filter'})

labels_x = repmat(" ", length(win_list) * 2, 1);
labels_x(2:2:end) = compose("%.2f", win_list' * 1e3);
set(gca,'XTickLabel', labels_x);

grid on
ylabel("Mean Absolute Error / Null Error")
xlabel("Window length [mm]")
title("Effect of window length on displacement estimation error")