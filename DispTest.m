addpath('lib/Auxiliar/')
addpath('lib/poly2D/')
addpath('lib/Field_II/')
rng(6942069)
field_init(0)
graf = 1;
tic()
%% Simulation parameters
% Siemens VF10-5 
% = [fs, c_c, f0, n_elem, elem_w, elv_focus, tx_focus, att, n_active] 
sim_data = [50e6, 1540, 6.67e6, 128, 0.3e-3, 20e-3, 10e-3, 70e-6, 64];

% Configurate simulation
simulation_config(sim_data);

% Configurate transducer 
trans_tx = create_transducer(sim_data);
trans_rx = create_transducer(sim_data);
sprintf("%.2f : Configuration ready", toc())

%% Place scatterers
% Scatterer density calculation
cell_vol = resolution_cell(sim_data, 0);
density = 5 / (cell_vol * 1e9); % [scatters/mm^3]

% Dispose scatterers at random
scat_pos = create_scatterers_3d([0, 6e-3], [0, 0.3e-3], [0, 6e-3], density);
%scat_pos = [3, 0, 3; 4, 0, 4] * 1e-3;
sprintf("%.2f : Scatterers ready", toc())

%% Generate pre-deformation sonogram
img_x = (0:0.05:8) * 1e-3;
[img_z, bmode_img] = bmode_image(trans_tx, trans_rx, sim_data, ...
                                         scat_pos, img_x, 1.2e3);
img_z = img_z - 0.62e-3; % Why?

sonograms = zeros(2, length(img_x), length(img_z));
sonograms(1, :, :) = bmode_img;
sprintf("%.2f : Pre B-Mode ready", toc())

%% Create dummy displacement
disp_x = (0:0.1:10) * 1e-3;
disp_y = (-1:1:1) * 1e-3;
disp_z = (0:0.1:15) * 1e-3;
[y_grid, x_grid, z_grid] = meshgrid(disp_y, disp_x, disp_z);

u_x = zeros(length(disp_x), 3, length(disp_z));
u_y = zeros(length(disp_x), 3, length(disp_z));
disp_mag = 0.1e-3;
%u_z = disp_mag * ones(length(disp_x), 3, length(disp_z));
u_z = disp_mag * exp(-((x_grid - 3e-3).^2 + (z_grid - 3e-3).^2) / (2 * (0.5e-3)^2));

%% Apply displacement to scatterers
new_scat_pos = apply_u(scat_pos, disp_x, disp_y, disp_z, u_x, u_y, u_z);
sprintf("%.2f : Displacement ready", toc())

%% Calculate post-deformation sonogram
[~, bmode_img] = bmode_image(trans_tx, trans_rx, sim_data, ...
                                 new_scat_pos, img_x, 1.2e3);
sonograms(2, :, :) = bmode_img;

% Terminate program and free memory 
xdc_free(trans_tx)
xdc_free(trans_rx)
field_end()
sprintf("%.2f : Post B-Mode ready", toc())

%% Estimate displacement
for win_len = [0.5, 1, 1.5, 2, 2.5] * 1e-3
    %win_len = 1.5e-3; % Window length [m]
    fine_prec = 0.1;   % Fine precission [%/%]
    [est_x, est_z, u_x_est, u_z_est] = estimate_u(img_x, img_z, sonograms,...
                                           win_len, 1.5e-4, fine_prec, 0);
    sprintf("%.2f : Disp. Estimation ready", toc())
    
    % Apply Filter to estimation
    u_est_filt = 3 * conv2(medfilt2(squeeze(u_z_est(1, :, :)), [7, 7]), ones(3, 3)/9, 'same');
    
    % Calculate error
    [u_z_err, total_error] = u_error(est_x, est_z, u_est_filt, ...
                                     disp_x, disp_z, squeeze(u_z(:, 2, :)));

    % Show results
    
    if ~graf
        fig = figure(1);
        fig.Position = [0, 0, 500, 500];
        imagesc(img_z, img_x, squeeze(sonograms(1, :, :)));
    %     hold on
    %     scatter(scat_pos(:, 3), scat_pos(:, 1), [],'red')
    %     hold off
        ylabel('Lateral Distance (X) [m]')
        xlabel('Depth (Z) [m]')
        title('First Sonogram')
        colormap(jet)
        colorbar
    
        fig = figure(2);
        fig.Position = [500, 0, 500, 500];
        imagesc(img_z, img_x, squeeze(sonograms(2, :, :)));
    %     hold on
    %     scatter(new_scat_pos(:, 3), new_scat_pos(:, 1), [],'red')
    %     hold off
        ylabel('Lateral Distance (X) [m]')
        xlabel('Depth (Z) [m]')
        title('Second Sonogram')
        colormap(jet)
        colorbar
    
        figure(3)
        scatter(scat_pos(:, 3), scat_pos(:, 1))
        hold on
        scatter(new_scat_pos(:, 3), new_scat_pos(:, 1))
        hold off
        axis ij
        axis([img_z(1), img_z(end), img_x(1), img_x(end)])
        title('Scatterer Position')
        xlabel('Depth (m)')
        ylabel('Lateral position (m)')
        legend({'Pre', 'Post'})
    end
    
    if graf
        fig = figure(4);
        fig.Position = [0, 0, 1500, 500];
        sgtitle('Z Displacement') 
    
        ax1 = subplot(1, 3, 1);
        fig_mesh1 = mesh(est_z, est_x, u_est_filt);
        zlim([0, 3 * disp_mag])
        clim([0, 3 * disp_mag])
        view(30, 30)
        title('Estimated')
    
        ax2 = subplot(1, 3, 2);
        fig_mesh2 = mesh(disp_z, disp_x, squeeze(u_z(:, 2, :)));
        title('Real')
    
        ax3 = subplot(1, 3, 3);
        fig_mesh3 = mesh(est_z, est_x, u_z_err);
        title(sprintf('Abs. Error: %.2f mm', total_error * 1e3))
        colorbar
    
        Link = linkprop([ax1, ax2, ax3],{'CameraUpVector', 'CameraPosition', 'CameraTarget', 'XLim', 'YLim', 'ZLim'});
        setappdata(gcf, 'StoreTheLink', Link);
        
        saveas(fig, sprintf("out_fig/est_d%.0f_win_%.1f.fig", density, win_len*1e3))
    
        for t_i = 2:size(u_z_est, 1)
            break
            pause(0.1)
            fig_mesh1.ZData = squeeze(u_z_est(t_i, :, :))';
            fig_mesh2.ZData = squeeze(u_z(t_i, :, :))';
        end
    end
end