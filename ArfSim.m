addpath('lib/Field_II/')
field_init(0)
clc

%% Configuramos simulación
fs = 100e6;  % Frecuencia de sampleo [Hz]
f0 = 6.67e6; % Frecuencia central de transductor [Hz]
set_sampling(fs);

% Configuramos propiedades acústicas
c_c = 1540;    % Velocidad del sonido [m/s]
rho = 1000;    % Densidad del medio [kg/m^3]
imp = rho * c_c; % Impedancia del medio

% Configuramos atenuación
alpha = 0.7 / (1e-2 * 1e6);             % Coeficiente de atenuación [dB/(m*Hz)]
absorcion = alpha * f0 / 20 * log(10);  % Coeficiente de absorción [Np/m]
set_field('att_f0', f0)                 % Frecuencia central
set_field('freq_att', alpha)            % 0.7 dB/ (cm * MHz)
set_field('att', f0 * alpha)            % Atenuación independiente de frec.
set_field('use_att', 1);

%% Configuración de Transductor (Siemens VF10-5)

n_elem = 128;   %  Número de elementos
width = 0.3e-3; % Ancho de elementos [m]
kerf = 0e-3;    % Distancia en x entre elementos [m]
height = 5e-3;  % Altura de elementos [m]

focus = [0 0 10] * 1e-3; % Foco electrónico inicial [m]

% Nr. elementos activos
n_act = ceil(1.3 * 10e-3 / (width+kerf));

%  Definimos el transductor
trans = xdc_focused_array(n_elem, width, height, kerf, 20e-3, ...
    1, 10, focus);

%figure(1)
%xdc_draw(trans)

% Creamos respuesta al impulso del transductor con una ventana Hanning
resp_imp = sin(2 * pi * f0 * (0:1/fs:4/f0));
resp_imp = resp_imp .* hanning(length(resp_imp))';
xdc_impulse(trans, resp_imp)

% Definimos voltaje de excitación del transductor
n_push = 300; % Número de ciclos en secuencia de push
volt_exc = sin(2 * pi * f0 * (0:1/fs:n_push/f0));
xdc_excitation(trans, volt_exc);

% Usamos baffle
xdc_baffle(trans, 1);

% figure(2)
% subplot(1, 2, 1)
% plot((0:length(resp_imp)-1)/fs, resp_imp)
% subplot(1, 2, 2)
% plot((0:length(volt_exc)-1)/fs, volt_exc)


%% Generamos grilla para calcular la intensidad acústica (0.1 mm)

dim_x = (-7:0.1:7) * 1e-3;
dim_z = (0:0.1:25) * 1e-3;

[X, Z] = ndgrid(dim_x, dim_z);
grilla = [X(:), zeros(length(X(:)), 1) , Z(:)];

%% Calculamos intensidad en cada punto

% Enfocamos transductor
xdc_center_focus(trans, [0, 0, 0]);
xdc_focus(trans, 0, focus);

% Activamos elementos
apodizacion = [zeros(1, floor((n_elem-n_act)/2)),...
               ones(1, n_act),...
               zeros(1, ceil((n_elem-n_act)/2))];
xdc_apodization(trans, 0, apodizacion);

% Calculamos campo de presión [Pa = N/m^2]
[presion, t_min] = calc_hp(trans, grilla);

%% Mostramos promedio temporal de intensidad [W/m^2 = 0.0001 W/cm^2]

Ispta = 1000 * 100^2; % Intensidad peak espacial, promedio temporal [W/m^2]
inten_prom = reshape(mean(presion.^2 / imp, 1), size(X))';
inten_prom = inten_prom / max(inten_prom, [], 'all') * Ispta;

% figure(3)
% mesh(X, Z, inten_prom')
% xlabel("Distancia lateral [m]")
% ylabel("Distancia axial [m]")
% title("Intensidad promedio-temporal [W/m^2]")
% 
% figure(4)
% imagesc(dim_x, dim_z, inten_prom)
% xlabel("Distancia lateral [m]")
% ylabel("Distancia axial [m]")
% title("Intensidad promedio-temporal [W/m^2]")
% axis image
% colorbar

%% Obtenemos fuerza donde la intensidad es alta (> 10% peak)

puntos_int = inten_prom > 0.1 * Ispta;
[ind_z, ind_x] = find(puntos_int);
f_rad_acus = (2 * absorcion * inten_prom / c_c) .* puntos_int;

%% Exportamos variables para FEniCSx
save('output/ARF.mat', 'f_rad_acus', 'ind_x', 'ind_z', 'dim_x', 'dim_z')

%% Mostramos resultados

load('output/ARF.mat', 'f_rad_acus', 'ind_x', 'ind_z', 'dim_x', 'dim_z')
[X, Z] = ndgrid(dim_x, dim_z);

figure(5)
mesh(X, Z, f_rad_acus')
xlabel("Distancia lateral [m]")
ylabel("Distancia axial [m]")
title("Fuerza de radiación acústica [N/m^3]")

%% Finalizamos
xdc_free(trans)
field_end
