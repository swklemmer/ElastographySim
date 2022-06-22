addpath('lib/Auxiliar/')
addpath('lib/poly2D/')
graf = 0;
E = 6000;       % Young's Modulus [Pa]
density = 28.7; % [scatters/cell]

win_len = 1.5e-3;  % Window size [m]
win_hop = 0.15e-3; % Hop between windows [m]
fine_prec = 0.1;   % Fine precission [%/%]

%% Estimate displacement from consecutive sonograms

% Simulation paramters
% Siemens VF10-5 
% = [fs, c_c, f0, n_elem, elem_w, elv_focus, tx_focus, att, n_active] 
sim_data = [50e6, 1540, 6.67e6, 128, 0.3e-3, 20e-3, 10e-3, 70e-6, 64];

% Import sonograms
load(sprintf('output/s_%d_d%.1f.mat', E, density),...
    'sonograms', 'img_x', 'img_z');
%%
if graf
    % Show first sonogram
    fig = figure(1);
    fig.Position = [0, 750, 750, 750];
    img = show_sonogram(img_x, img_z, sonograms,...
        0, [10, 15]*1e-3, 1, 'First Sonogram');
    show_sequence(sonograms, img)
end

%% Estimate displacement
[est_x, est_z, u_x_est, u_z_est] = estimate_u(img_x, img_z, ...
                                        sonograms, win_len, win_hop, ...
                                        fine_prec, graf);

%% Save results
save(sprintf('output/u_est_%d_d%d.mat', E, density),...
    'u_x_est', 'u_z_est', 'est_x', 'est_z')

%% Import real and estimated displacement
[t_dim, x_dim, y_dim, z_dim, ~, u_z] = ...
    load_mid_u(sprintf("input/u_%d.h5", E));

load(sprintf('output/u_est_%d_d%.1f.mat', E, density),...
    'u_x_est', 'u_z_est', 'est_x', 'est_z')

%% Apply filtering
filt_sz = 5;
u_filt = u_z_est;
for t = 1:size(u_filt, 1)
    u_filt(t, :, :) = conv2(squeeze(u_z_est(t, :, :)),...
                            ones(filt_sz)/filt_sz^2, 'same');
end

%% Normalize displacements
u_z = u_z / max(u_z, [], 'all');
u_filt = u_filt / max(u_filt, [], 'all');

%% Show results

if graf
    [u_err, mean_err] = estimation_error(...                               
                                x_dim, z_dim, squeeze(u_z(2, :, :)),...
                                est_x, est_z, squeeze(u_filt(1, :, :)), ...
                                [8, 15]*1e-3);
    fig = figure(1);
    fig.Position = [0, 0, 1000, 300];
    [img_real, img_est, img_err] = show_disp_error(...
                                x_dim, z_dim, squeeze(u_z(2, :, :)), ...
                                est_x, est_z, squeeze(u_filt(1, :, :)), ...
                                u_err, 0, [8, 15]*1e-3);
 
    for t = 2:size(u_filt, 1)
        pause(0.1)
        [u_err, ~] = estimation_error(...                               
                                x_dim, z_dim, squeeze(u_z(t+1, :, :)),...
                                est_x, est_z, squeeze(u_filt(t, :, :)), ...
                                [8, 15]*1e-3);

        img_real.ZData = squeeze(u_z(t + 1, :, :));
        img_est.ZData  = squeeze(u_filt(t, :, :));
        img_err.ZData  = u_err;
    end
end
