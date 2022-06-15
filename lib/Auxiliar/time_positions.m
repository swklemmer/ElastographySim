function img_z = time_positions(sim_data, z_crop)
%TIME_POSITIONS Calculate z positions for sonograms

fs = sim_data(1);        % Sampling Frequency [Hz]
c_c = sim_data(2);       % Sound Speed [m/s]

img_z = (1:z_crop) * c_c / fs / 2;

end
