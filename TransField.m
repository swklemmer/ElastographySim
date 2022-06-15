addpath('lib/Field_II/')
addpath('lib/Auxiliar/')
rng(6942069)
field_init(0)

% Simulation paramters

% Siemens VF10-5 = [fs, c_c, f0, n_elem, elem_w, elv_focus, tx_focus] 
sim_data = [50e6, 1540, 6.67e6, 128, 0.3e-3, 10e-3, 10e-3];

% Recover simulation data
fs = sim_data(1);        % Sampling Frequency [Hz]
c_c = sim_data(2);       % Sound Speed [m/s]
n_elem = sim_data(4);    % Number of elements
elem_w = sim_data(5);    % Element width [m]
tx_focus = sim_data(7);  % Transmit focus depth [m]
n_active = 63;

% Read info from h5 file
[~, x_dim, y_dim, z_dim, ~, ~, ~] = load_u("input/u_4000.h5");

% Configurate simulation
simulation_config(sim_data);

% Transducer configuration
trans_tx = create_transducer(sim_data);
trans_rx = create_transducer(sim_data);

% Find center element
center_elem = n_elem / 2;

% Define active elements around it
left_elem = center_elem - (n_active - 1) / 2;
right_elem = center_elem + (n_active - 1) / 2;
active_elem = max(1, left_elem):min(right_elem, n_elem);

% Activate elements through apodization
hann_win = hanning(n_active);
apodization = zeros(1, n_elem);

if left_elem <= 0
    apodization(active_elem) = hann_win(2-left_elem:end);
elseif right_elem >= n_elem
    apodization(active_elem) = hann_win(1:end-(right_elem-n_elem));
else
    apodization(active_elem) = hann_win;
end    

xdc_apodization(trans_tx, 0, apodization);
xdc_apodization(trans_rx, 0, apodization);

% Focus transducer
xdc_center_focus(trans_tx, [0, 0, 0]);
xdc_focus(trans_tx, 0, [0, 0, tx_focus]);

xdc_center_focus(trans_rx, [0, 0, 0]);
xdc_dynamic_focus(trans_rx, 0, 0, 0);

% Calculate pulse-echo field
pul_ech = zeros(length(x_dim), length(y_dim), length(z_dim));

for z = 1:length(z_dim) 
    for y = 1:length(y_dim)
        for x = 1:length(x_dim)
            [hhp, start_time] = calc_hhp(trans_tx, trans_rx, [x_dim(x), y_dim(y), z_dim(z)]);
            pul_ech(x, y, z) = max(hhp.^2, [], 'all');
        end
    end
end

pul_ech = pul_ech / max(pul_ech, [], 'all');

%% Show field
fig_mesh = mesh(x_dim, y_dim, squeeze(pul_ech(:, :, 1))');
view(80, 5)
xlabel('Lateral distance (m)')
ylabel('Elevation distance (m)')
zlabel('Acoustic Pressure (U.A)')
zlim([0, 1])
clim([0, 1])
colorbar

for z = 2:length(z_dim)
    pause(0.1)
    fig_mesh.ZData = squeeze(pul_ech(:, :, z))';
    title(sprintf("Normalized pressure at z = %.1f mm", 1e3 * z_dim(z)));
end

%% Terminate program and free memory 
xdc_free(trans_tx)
xdc_free(trans_rx)
field_end