function cell_vol = resolution_cell(sim_data, graf)
%RESOLUTION_CELL Calculate de -6dB elliptic cylinder volume of a given
%transducer configuration

% Recover simulation data
fs = sim_data(1);        % Sampling Frequency [Hz]
c_c = sim_data(2);       % Sound Speed [m/s]
f0 = sim_data(3);        % Transducer center frequency [Hz]
n_elem = sim_data(4);    % Number of elements
tx_focus = sim_data(7);  % Transmit focus depth [m]
alpha = sim_data(8);     % Attenuation coefficient [dB/m/Hz]
n_active = sim_data(9);  % Nr. of active elements

% Create analysis volume
x_dim = (-2:0.02:2)*1e-3;
y_dim = (-2:0.02:2)*1e-3;
z_dim = (-2:0.02:2)*1e-3 + tx_focus;

% Configurate simulation
simulation_config(sim_data);

% Transducer configuration
trans_tx = create_transducer(sim_data);
trans_rx = create_transducer(sim_data);

% Find center element
center_elem = n_elem / 2;

% Define active elements around it
left_elem = center_elem - n_active / 2 + 1;
right_elem = center_elem + n_active / 2;
active_elem = max(1, left_elem):min(right_elem, n_elem);

% Activate receive elements through apodization
hann_win = hanning(n_active);
apodization = zeros(1, n_elem);

if left_elem <= 0
    apodization(active_elem) = hann_win(2-left_elem:end);
elseif right_elem >= n_elem
    apodization(active_elem) = hann_win(1:end-(right_elem-n_elem));
else
    apodization(active_elem) = hann_win;
end    

xdc_apodization(trans_rx, 0, apodization);
xdc_apodization(trans_tx, 0, apodization);

% Focus transducer
xdc_center_focus(trans_tx, [0, 0, 0]);
xdc_focus(trans_tx, 0, [0, 0, tx_focus]);

xdc_center_focus(trans_rx, [0, 0, 0]);
xdc_focus(trans_rx, 0, [0, 0, tx_focus]);

% Calculate pulse-echo field in x direction
pul_ech_x = zeros(length(x_dim), 1);

for x = 1:length(x_dim)
    [hhp, start_time] = calc_hhp(trans_tx, trans_rx, [x_dim(x), 0, tx_focus]);
    pul_ech_x(x) = max(abs(hilbert(hhp.^2)), [], 'all');
end

% Calculate pulse-echo field in y direction
pul_ech_y = zeros(length(y_dim), 1);

for y = 1:length(y_dim)
    [hhp, ~] = calc_hhp(trans_tx, trans_rx, [0, y_dim(y), tx_focus]);
    pul_ech_y(y) = max(abs(hilbert(hhp.^2)), [], 'all');
end

% Calculate pulse-echo field in z direction
pul_ech_z = zeros(length(z_dim), 1);

for z = 1:length(z_dim)
    [hhp, ~] = calc_hhp(trans_tx, trans_rx, [0, 0, z_dim(z)]);

    % Apply Time-Gain-Compensation
    hhp_zp = [zeros(round(start_time * fs), 1); hhp];
    att_dB = alpha * f0 * (1:length(hhp_zp)) * c_c / fs / 2;
    hhp_zp = hhp_zp' .* 10.^(att_dB / 20);
    pul_ech_z(z) = max(abs(hilbert(hhp_zp.^2)), [], 'all');
end

% Normalize every distribution
pul_ech_x = pul_ech_x / max(pul_ech_x, [], 'all');
pul_ech_y = pul_ech_y / max(pul_ech_y, [], 'all');
pul_ech_z = pul_ech_z / max(pul_ech_z, [], 'all');

% Calculate volume
x_0 = find(pul_ech_x > 10^(-6/20), 1, "first");
x_f = find(pul_ech_x > 10^(-6/20), 1, "last");
y_0 = find(pul_ech_y > 10^(-6/20), 1, "first");
y_f = find(pul_ech_y > 10^(-6/20), 1, "last");
z_0 = find(pul_ech_z > 10^(-6/20), 1, "first");
z_f = find(pul_ech_z > 10^(-6/20), 1, "last");

cell_vol = pi / 4 * (x_dim(x_f) - x_dim(x_0)) * ...
                    (y_dim(y_f) - y_dim(y_0)) * ...
                    (z_dim(z_f) - z_dim(z_0));

% Free memomry
xdc_free(trans_tx);
xdc_free(trans_rx);

% Show results
if graf
fig = figure(1);
plot(x_dim * 1e3, pul_ech_x);
fig.Position = [0, 200, 300, 200];
yline(10^(-6/20), '--')
xlabel("Lateral distance [mm]")
ylabel("Normalized transducer response")
title("-6dB resolution cell in X direction")
grid on

fig = figure(2);
plot(y_dim * 1e3, pul_ech_y);
fig.Position = [400, 200, 300, 200];
yline(10^(-6/20), '--')
xlabel("Elevation distance [mm]")
ylabel("Normalized transducer response")
title("-6dB resolution cell in Y direction")
grid on

fig = figure(3);
plot(z_dim * 1e3, pul_ech_z);
fig.Position = [800, 200, 300, 200];
yline(10^(-6/20), '--')
xlabel("Axial distance [mm]")
ylabel("Normalized transducer response")
title("-6dB resolution cell in Z direction")
grid on
end
end

