addpath('lib/Field_II/')
addpath('lib/Auxiliar/')
rng(6942069)
field_init(0)

% Simulation paramters
% Siemens VF10-5 
% = [fs, c_c, f0, n_elem, elem_w, elv_focus, tx_focus, att, n_active] 
sim_data = [50e6, 1540, 6.67e6, 128, 0.3e-3, 20e-3, 10e-3, 70e-6, 64];

% Configurate simulation
simulation_config(sim_data);

% Transducer configuration
trans_tx = create_transducer(sim_data);
trans_rx = create_transducer(sim_data);

% Import FEM model info
E = 4000; % [Pa]
fem_file = sprintf("input/u_%d.h5", E);

% Read info from h5 file
[t_dim, x_dim, y_dim, z_dim, u_x, u_y, u_z] = load_u(fem_file);

% Allocate B-mode adquisition data
d_list = [1, 2, 5, 10] * 1e-1; % [scat/cell]

z_crop = 1.5e3;
img_x = x_dim(1):0.1e-3:x_dim(end); % 8 divisions per res. cell
sonograms = zeros(length(d_list), length(img_x), z_crop);

for d = 1:length(d_list)
    tic()
    % Scatterer density calculation
    cell_vol = resolution_cell(sim_data, 0);
    density = d / (cell_vol * 1e9); % [scatters/mm^3]

    % Dispose scatterers at random
    scat_pos = create_scatterers_3d(x_dim, y_dim(1:4), z_dim, density);
    
    % Generate image
    [img_z, bmode_img] = bmode_image(trans_tx, trans_rx, sim_data, ...
                                     scat_pos, img_x, 2.5e3);

    % Save image
    sonograms(d, :, :) = bmode_img(:, 1:z_crop);

    it_time = toc();
    fprintf("\nProgreso = %.1f %%\nRestante = %.1f min\n", ...
                                        d / length(d_list) * 100,...
                                    it_time * (length(d_list) - d) / 60)
end

% Save results
img_z = img_z(1:z_crop);
img_z = img_z - 0.62e-3; % WHY?
save('out_dens/s_density.mat', 'sonograms', 'img_x', 'img_z', 'd_list');

% Terminate program and free memory 
xdc_free(trans_tx)
xdc_free(trans_rx)
field_end()

%% Load Sonograms
load('out_dens/s_density.mat', 'sonograms', 'img_x', 'img_z');

%% Show intensity histograms

% Create Rayleigh pdf
std_ = linspace(0, 6, 60);
p_rayleigh = std_ .* exp( - std_.^2 / 2);

% Plot all histograms together
fig = figure(1);
fig.Position = [500, 500, 500, 250];
hold on
for d = 1:length(d_list)
    phantom_crop = squeeze(sonograms(d, :, 2e-3 < img_z & img_z < 17e-3));
    inten_sd = std(phantom_crop, 0, 'all');
    histogram(phantom_crop / inten_sd, 'Normalization', 'pdf')
end
plot(std_, p_rayleigh, 'LineWidth', 3, 'LineStyle', '--')
hold off
xlim([-0.1, 6])
xlabel("Amplitude / STD")
ylabel("Amplitude PDF")
title("Effect of scatterer density on B-Mode Histograms")
legend(['1 scat/mm^3', '2 scat/mm^3', '5 scat/mm^3',...
    '10 scat/mm^3', "Rayleigh distribution"])
grid on
