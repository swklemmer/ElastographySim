function [v_est, v_theo, v_error] = ...
    TTP_estimate_v(est_t, est_x, est_z, est_u, E, sim_data, x_span, graf)
%TTP_ESTIMATE_V Estimates shear wave speed [m/s] for a given sonogram
%sequencue s, Young's Modulus E [kPa]
G = E/3; % Pa 
rho = 1000;
v_theo = sqrt(G/rho);
z_focus = sim_data(7); % ARF focal depth [m]

% Average over focal depth
delta_z = 1e-3; % Window width [mm]
z_start = find(est_z >= z_focus - delta_z, 1, 'first');
z_end = find(est_z >= z_focus + delta_z, 1, 'first');

u_focal = sum(est_u(:, :, z_start:z_end) / (z_end-z_start), 3);

% Define lateral ROI
min_x = find(est_x >= x_span(1), 1, 'first');
max_x = find(est_x >= x_span(2), 1, 'first');

% Shear Wave Speed (Time-To-Peak Reconstruction)
ttpeak = zeros(1, length(est_x));
hres_dim = linspace(est_t(1), est_t(end), 2^12);

for x = 1:length(est_x)
    % Interpolate data
    hres_data = interp1(est_t(1:size(u_focal, 1)), u_focal(:, x), hres_dim, 'cubic');

    % Find temporal position of peak
    [~, peak_ind] = max(hres_data);
    ttpeak(x) = hres_dim(peak_ind);
end

% Fit linear regression
fitted_line = polyfit(est_x(min_x:max_x), ttpeak(min_x:max_x), 1);
v_est = 1/fitted_line(1);

v_error = RSE(v_theo, v_est);

% Show results
if graf
    fprintf("Youngs = %d\nTh. Ct = %.3f\nFEM Ct = %.3f\nError  = %.1f %%", ...
    E, v_theo, v_est, v_error)

    % Plot TTP
    f = figure();
    f.Position = [0, 0, 300, 300];
    plot(est_x(min_x:max_x), ttpeak(min_x:max_x))
    hold on
    plot(est_x(min_x:max_x), est_x(min_x:max_x) / v_theo, '--')
    hold off
    grid on
    xlabel("Lateral position [m]")
    ylabel("Time to Peak [s]")
    title("Time to Peak Reconstruction")
    legend(["FEM", "Theorical"])

    % Plot Peak Profile for different times
    f = figure();
    f.Position = [0, 300, 300, 300];
    plot(est_t(1:size(u_focal, 1)), u_focal(:, 1) / max(u_focal(:, 1), [], 'all'))
    hold on
    for x = min_x:2:max_x
        plot(est_t(1:size(u_focal, 1)), u_focal(:, x) / max(u_focal(:, x), [], 'all'))
    end
    hold off
    ylim([-0.1, 1.1])
    xlabel('Time (s)')
    ylabel('Norm. Z-Displacement')
    grid
    title('Peak Profile through Time')
end
end
