addpath('lib/Auxiliar/')
addpath('lib/poly2D/')
addpath('lib/Field_II/')
rng(6942069)
graf = 0;
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
density = 6 / (cell_vol * 1e9); % [scatterers/mm^3]

% Iterate over various estimation window lengths
x_span = [1.3, 2.7] * 1e-3;
span_list = [2.7] * 1e-3;

% Estimate shear wave speed for various files
E_list = [4000, 6000];

errors = zeros(length(span_list), length(E_list));

for e = 1:length(E_list)
    E = E_list(e); % Young's Mod. [Pa]     

   % Load estimated displacement
    load(sprintf("output/u_est_%d_d%.1f.mat", E, density), 'est_x', 'est_z', 'u_z_est');
    [est_t, ~, ~, ~, ~, ~] = load_mid_u(sprintf("input/u_%d.h5", E));

    for i = 1:length(span_list)

        x_span(2) = span_list(i);
    
        % Estimate shear wave speed
        [v_est, v_theo, v_error] = ...
                              TTP_estimate_v(est_t, est_x, est_z, u_z_est, ...
                                             E, sim_data, x_span, graf);
    
        fprintf("\nE: %d, xspan: %.2g, %.2g, error = %.1f %%\n", ...
            E, x_span(1), x_span(2), v_error)

        errors(i, e) = v_error;
    end
end

%% Show results

% fig = figure(1);
% plot(span_list, errors(:, 1))
% hold on
% plot(span_list, errors(:, 2))
% hold off