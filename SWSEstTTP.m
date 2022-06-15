%% Load aux functions
addpath("lib/Auxiliar/")
close all
graf = 1;

for E = [4000, 6000, 8500, 10000]

% FEM parameters
G = E/3; % Pa 
rho = 1000;
C_theo = sqrt(G/rho);

% Import FEM data
fem_file = sprintf("input/u_%d.h5", E);
[t_dim, x_dim, ~, z_dim, ~, u_mat_z] = load_mid_u(fem_file);

%% Average over focal depth

% Focus
z_center = 10e-3; % m
delta_z = 1e-3; %mm

z_start = find(z_dim >= z_center - delta_z, 1, 'first');
z_end = find(z_dim >= z_center + delta_z, 1, 'first');

u_focal = sum(u_mat_z(:, :, z_start:z_end) / (z_end-z_start), 3);

% Define lateral ROI
min_x = find(x_dim >= 0.5e-3, 1, 'first');
max_x = find(x_dim >= 3.5e-3, 1, 'first');


%% Shear Wave Speed (Time-To-Peak Reconstruction)

ttpeak = zeros(1, length(x_dim));
hres_dim = linspace(t_dim(1), t_dim(end), 2^12);

for x = 1:length(x_dim)
    % Interpolate data
    hres_data = interp1(t_dim, u_focal(:, x), hres_dim, 'cubic');

    % Find temporal position of peak
    [~, peak_ind] = max(hres_data);
    ttpeak(x) = hres_dim(peak_ind);
end



% Fit linear regression
fitted_line = polyfit(x_dim(min_x:max_x), ttpeak(min_x:max_x), 1);
C_fem = 1/fitted_line(1);

sprintf("Youngs = %d\nTh. Ct = %.3f\nFEM Ct = %.3f\nError  = %.1f %%", ...
    E, C_theo, C_fem, RSE(C_theo, C_fem))

if graf
    
    % Plot TTP
    f = figure(1);
    f.Position = [0, 0, 300, 300];
    plot(x_dim(min_x:max_x), ttpeak(min_x:max_x))
    hold on
    plot(x_dim(min_x:max_x), x_dim(min_x:max_x) / C_theo, '--')
    hold off
    grid on
    xlabel("Lateral position [m]")
    ylabel("Time to Peak [s]")
    title("Time to Peak Reconstruction")
    legend(["FEM", "Theorical"])

    % Plot Peak Profile for different times
    f = figure(2);
    f.Position = [0, 300, 300, 300];
    plot(t_dim, u_focal(:, 1))
    hold on
    for x = 2:1:max_x
        plot(t_dim, u_focal(:, x))
    end
    hold off
    xlabel('Time (s)')
    ylabel('Z-Displacement (m)')
    grid
    title('Peak Profile through Time')
end

%% Shear Wave Speed (Time-To-Peak Slope Reconstruction)

% Particle Velocity
v_focal = (u_focal - [zeros(1, length(x_dim)); u_focal(1:end-1, :)]) /...
    (t_dim(2));

ttpeak_s = zeros(1, length(x_dim));
hres_dim = linspace(t_dim(1), t_dim(end), 2^12);

for x = 1:length(x_dim)
    % Interpolate data
    hres_data = interp1(t_dim, v_focal(:, x), hres_dim, 'cubic');

    % Find temporal position of peak
    [~, peak_ind] = max(hres_data);
    ttpeak_s(x) = hres_dim(peak_ind);
end

if graf
    
    % Plot TTPS
    f = figure(3);
    f.Position = [300, 0, 300, 300];
    plot(x_dim(min_x:max_x), ttpeak_s(min_x:max_x))
    hold on
    plot(x_dim(min_x:max_x), x_dim(min_x:max_x) / sqrt(G/rho), '--')
    hold off
    grid on
    xlabel("Lateral position [m]")
    ylabel("Time to Peak [s]")
    title("Time to Peak Slope Reconstruction")
    legend(["FEM", "Theorical"])

    % Plot Peak Profile for different times
    f = figure(4);
    f.Position = [300, 300, 300, 300];
    plot(t_dim, v_focal(:, 1))
    hold on
    for x = 2:1:max_x
        plot(t_dim, v_focal(:, x))
    end
    hold off
    xlabel('Time (s)')
    ylabel('Particle Velocity (m/s)')
    grid
    title('Velocity Profile through Time')
end

waitforbuttonpress();
end

