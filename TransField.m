addpath('lib/Field_II/')
addpath('lib/Auxiliar/')
rng(6942069)
field_init(0)
graf = 1;

%% Simulation paramters
% Siemens VF10-5 = [fs, c_c, f0, n_elem, elem_w, elv_focus, tx_focus, att] 
sim_data = [50e6, 1540, 6.67e6, 128, 0.3e-3, 20e-3, 10e-3, 70e-6, 64];

%% Calculate resolution cell volume
cell_vol = resolution_cell(sim_data, graf);
fprintf("\n\nRes. Cell. Volume = %.2f mm^3\n\n", cell_vol * 1e9);

%% Terminate program
field_end