%% Load aux functions
addpath("lib/Auxiliar/")
addpath("Deprecated/EqBasedReconstruction/")

youngs = 10000;
fem_file = sprintf("input/u_%d.h5", youngs);
graf = 1;

%% Importamos una respuesta real

[t_dim, x_dim, y_dim, z_dim, u_matrix] = load_mid_u(fem_file);

% Parámetros reales
G_1 = 10; % Pa
G_inf = youngs/3; % Pa 
beta = 20000; % Hz

%% Reconstrucción por 2D-FT

% Truncation frecuencies
f_t_tr = 1400; % Hz
f_s_tr = 1000; % 1/m

% Distancia del foco
z_central = 10e-3; % m
delta_z = 0.5e-3; %mm

z_start = find(z_dim >= z_central - delta_z, 1, 'first');
z_end = find(z_dim >= z_central + delta_z, 1, 'first');

% Padding de FFT
pad_s = 2^12;
pad_t = 2^12;

%% Velocidad de partícula
% Sumamos desplazamiento en un segmento para aumentar SNR
esp_tiempo = diff(squeeze(sum(u_matrix(:, :, z_start:z_end), 3)), 1, 1)...
    ./ diff(t_dim)';
esp_k = abs(flip(fft2(esp_tiempo', pad_s, pad_t), 1));

% Región de interés
ds = x_dim(2) - x_dim(1);
frec_s = linspace(0, 1/(2 * ds), pad_s);
fs_max = find(frec_s > f_s_tr, 1, 'first');

dt = t_dim(2) - t_dim(1);
frec_t = linspace(0, 1/(2 * dt), pad_t);
ft_max = find(frec_t > f_t_tr, 1, 'first');

% Truncamos
esp_k_trunc = esp_k(1:fs_max, 1:ft_max);

% Máximo por frecuencia
[~, max_k] = max(esp_k_trunc, [], 1);

if graf
    % Graficamos máxima amplitud en frecuencia espacial-temporal
    figure(1);
    imagesc(frec_t(1:ft_max), frec_s(1:fs_max),  esp_k_trunc)
    hold on
    plot(frec_t(1:ft_max), frec_s(max_k), 'ro')
    hold off
    set(gca, 'YDir', 'normal')
    xlabel('Frecuencia (Hz)')
    ylabel('Frecuencia espacial (1/m)')
    title('Espacio k')
end

%% Graficamos Estimación de velocidad versus frecuencia

% Inferimos velocidad cortante de FEM
cs_fem = frec_t(1:ft_max) ./ frec_s(max_k);

% Filtramos datos extraños
%[~, frec_est_min] = min(cs_est);

%frec_est = frec_t(frec_est_min:floor(ft_max/3)) * 3;
%cs_est = cs_est(frec_est_min:floor(ft_max/3));

% Inferimos módulo de Young de FEM
G_est = G_cs(cs_fem);

% Respuesta mecánica teórica
[f_teo, G_teo, cs_teo] = zener_model(f_t_tr, G_1, G_inf, beta);

if graf
    figure(2);
    yyaxis left
    plot(f_teo, cs_teo, '--');
    hold on
    plot(frec_t(1:ft_max), cs_fem, '-');
    hold off
    grid on
    ylabel('Velocidad cortante [m/s]')
    ylim([0, 5])
    
    yyaxis right
    plot(f_teo, abs(G_teo)*1e-3, '--');
    hold on 
    plot(frec_t(1:ft_max), G_est*1e-3, '-');
    hold off
    grid on
    xlabel('Frecuencia [Hz]')
    xlim([0, f_t_tr])
    ylabel('|Módulo cortante| [kPa]')
    ylim([-4, 6])
    
    title('Reconstrucción de propiedades mecánicas')
    legend({"$c_s$ te\'orico", '$c_s$ FEM', "G te\'orico", 'G FEM'},...
        'position', [.7, .2, 0.1, 0.2], 'interpreter', 'Latex')
end

%% Encontramos parámetros mecánicos
% x1 = G_1, x2 = G_inf, x3 = beta

% param_est = lsqcurvefit(@(x, xdata) cs_G(x, xdata),...
%     [G_1, G_inf, beta], frec_est, double(cs_est),...
%     [5e3, 1e3, 1e3], [25e3, 10e3, 10e3]);
% 
% fprintf(['          Real  |  Estimado  |  Error\n',...
%         'Ginf: %8.1f  |  %8.1f  |  %2.2f%%\n',...
%         'eta : %8.2f  |  %8.2f  |  %2.2f%%\n'],...
%     G_inf, param_est(2), RSE(G_inf, param_est(2)), ...
%     G_1/beta, param_est(1)/param_est(3),...
%     RSE(G_1/beta, param_est(1)/param_est(3)))
% 
% 
% f_recon = linspace(0, f_tr, 100);
% 
% if graf
%     figure(3)
%     plot(f_recon, cs_G(param_est, f_recon))
%     hold on
%     plot(frec_est, cs_est, '-');
%     plot(f_real, cs_real, '--');
%     hold off
%     grid on
%     legend('Reconstruido', 'Observado', 'Real')
%     title('|Módulo cortante|')
% end
