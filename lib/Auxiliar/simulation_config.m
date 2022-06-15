function simulation_config(sim_data)
%SIMULATION_CONFIG Sets up Field II with following parameters
% Attenuation: 0.7 dB/ (cm*Hz)

% Returns: Acoustic Impedance and Absorption
fs = sim_data(1);  % Sampling frequency [Hz]
set_sampling(fs);

f0 = sim_data(3);  % Transducer center frequency [Hz]

% Configurate Attenuation
alpha = sim_data(8);                    % Attenuation Coeficient [dB/(m*Hz)]
%absorption = alpha * f0 / 20 * log(10);% Absorption Coeficient [Np/m]
set_field('att_f0', f0)                 % Absorption center frecuency
set_field('freq_att', alpha)            % 0.7 dB/ (cm * MHz)
set_field('att', f0 * alpha)            % Frequency ndependent attenuation
set_field('use_att', 1);

end

