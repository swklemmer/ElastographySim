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

% Scatterer density calculation
cell_vol = resolution_cell(sim_data, 0);
density = 4 / (cell_vol * 1e9); % [scatters/mm^3]

%%
for E = [4000, 6000, 8500, 10000]

    % Import FEM model info
    fem_file = sprintf("input/u_%d.h5", E);
    
    % Read info from h5 file
    [t_dim, x_dim, y_dim, z_dim, u_x, u_y, u_z] = load_u(fem_file);

    % Dispose scatterers at random
    scat_pos = create_scatterers_3d(x_dim, y_dim(1:4), z_dim, density);
    %scat_pos = [(0:0.8:5)', zeros(7, 1), ones(7, 1) * 10] * 1e-3;

    % Simulate displacement and B-mode adquisition
    z_crop = 1.5e3;
    img_x = x_dim(1):0.1e-3:x_dim(end); % 8 divisions per res. cell
    sonograms = zeros(length(t_dim), length(img_x), z_crop);
    
    for t = 1:length(t_dim)
        tic()
    
        % Apply displacement to scatterers
        new_scat_pos = apply_u(scat_pos, x_dim, y_dim, z_dim, ...
                                            squeeze(u_x(t, :, :, :)), ...
                                            squeeze(u_y(t, :, :, :)), ...
                                            squeeze(u_z(t, :, :, :)));
    
        % Generate image
        [img_z, bmode_img] = bmode_image(trans_tx, trans_rx, sim_data, ...
                                         new_scat_pos, img_x, 2.5e3);
    
        % Save image
        sonograms(t, :, :) = bmode_img(:, 1:z_crop);
    
        it_time = toc();
        fprintf("\nProgreso = %.1f %%\nRestante = %.1f min\n", ...
                                            t / length(t_dim) * 100,...
                                        it_time * (length(t_dim) - t) / 60)
    end
    % Save results
    img_z = img_z(1:z_crop);
    img_z = img_z - 0.62e-3; % WHY?
    save(sprintf('output/s_%d_d%.1f.mat', E, density),...
        'sonograms', 'img_x', 'img_z');
end

% Terminate program and free memory 
xdc_free(trans_tx)
xdc_free(trans_rx)
field_end()

%% Load Sonograms
E = 10000;    % Young's Modulus [Pa]
density = 4 / (cell_vol * 1e9); % [scatters/mm^3]
load(sprintf('output/s_%d_d%.1f.mat', E, density), 'sonograms', 'img_x', 'img_z');

%% Apply logarithmic compression
gamma = 1e23;
sonograms_log = log(1 + gamma * sonograms);

% Show results
fig = figure(1);
fig.Position = [0, 500, 500, 250];
img = show_sonogram(img_x, img_z, sonograms_log,...
    0, [10, 15]*1e-3, 1, 'Sonogram Sequence');

show_sequence(sonograms_log, img)
