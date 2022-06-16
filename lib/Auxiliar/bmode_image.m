function [img_z, bmode_img] = bmode_image(trans_tx, trans_rx, sim_data, ...
                                          scat_pos, img_x, t_end)
%BMODE_IMAGE Simulate B-Mode image for a given position of A-lines. It
%automatically finds the active elements, their apodization and focus.

% Recover simulation data
fs = sim_data(1);        % Sampling Frequency [Hz]
v_c = sim_data(2);       % Sound Speed [m/s]
f0 = sim_data(3);        % Transducer center frequency [Hz]
n_elem = sim_data(4);    % Number of elements
elem_w = sim_data(5);    % Element width [m]
tx_focus = sim_data(7);  % Transmit focus depth [m]
alpha = sim_data(8);     % Attenuation coefficient [dB/m/Hz]
n_active = sim_data(9);  % Active elements during reception

% Prealocate data
rf_data = zeros(length(img_x), t_end);
hann_win = hanning(n_active);

for i = 1:length(img_x)

    % Find center element
    center_elem = n_elem / 2 + floor(img_x(i) / elem_w + 0.5);

    % Define active elements around it
    left_elem = center_elem - n_active / 2 + 1;
    right_elem = center_elem + n_active / 2;
    active_elem = max(1, left_elem):min(right_elem, n_elem);

    % Activate receive elements through apodization
    apodization = zeros(1, n_elem);
    
    if left_elem <= 0
        apodization(active_elem) = hann_win(2-left_elem:end);
    elseif right_elem >= n_elem
        apodization(active_elem) = hann_win(1:end-(right_elem-n_elem));
    else
        apodization(active_elem) = hann_win;
    end    

    xdc_apodization(trans_rx, 0, apodization);

    % Activate single transmit element
    apodization = zeros(1, n_elem);
    apodization(center_elem) = 1;
    xdc_apodization(trans_tx, 0, apodization);

    % Focus transducer
    xdc_center_focus(trans_tx, [img_x(i), 0, 0]);
    xdc_focus(trans_tx, 0, [img_x(i), 0, tx_focus]);

    xdc_center_focus(trans_rx, [img_x(i), 0, 0]);
    xdc_dynamic_focus(trans_rx, 0, 0, 0);
    
    % Simulate echos
    [v_signal, t_start] = calc_scat(trans_tx, trans_rx, scat_pos, ...
                                    ones(size(scat_pos, 1), 1));

    % Zero-pad beginnning to account for 't_start'
    v_zerop = [zeros(round(t_start * fs), 1) ; v_signal];

    % Save read-out
    rf_data(i, 1:length(v_zerop)) = v_zerop;
end

% Calculate depth dimension
img_z = (1:t_end) * v_c / fs / 2;

% Calculate depth dependent attenuation  
att_dB = alpha * f0 * img_z;

% Demodulate signal
bmode_img = zeros(length(img_x), t_end);

for i = 1:length(img_x)
    % Take Hilbert Transform of every line
    bmode_img(i, :) = abs(hilbert(squeeze(rf_data(i, :))));

    % Apply Time-Gain-Compensation
    bmode_img(i, :) = bmode_img(i, :) .* 10.^(att_dB / 20);
end
end