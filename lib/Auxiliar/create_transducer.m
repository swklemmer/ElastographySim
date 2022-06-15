function trans = create_transducer(sim_data)
%STANDARD_TRANSDUCER Creates a Field II Transducer object based on the 
% Siemens VF10-5 Model.

n_elem = sim_data(4);  % Number of elements
width = sim_data(5);   % Element width [m]
kerf = 0e-3;           % Distance between elements [m]
height = 5e-3;         % Element height [m]

focus = [0 0 10] * 1e-3; % Initial electronic focus [m]
elv_focus = sim_data(6); % Elevation Focus [m]

% Define transducer object
trans = xdc_focused_array(n_elem, width, height, kerf, elv_focus, ...
                          1, 10, focus);

% Sampling and central frecuency
fs = sim_data(1);
f0 = sim_data(3);

% Create transducer's impulse response using a Hanning window
imp_res = sin(2 * pi * f0 * (0:1/fs:2/f0));
imp_res = imp_res .* hanning(length(imp_res))';
xdc_impulse(trans, imp_res)

% Define transducer's input voltaje signal
pulser_v = sin(2 * pi * f0 * (0:1/fs:2/f0));
xdc_excitation(trans, pulser_v);

% Activate baffle
xdc_baffle(trans, 1);
end

