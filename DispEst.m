addpath('lib/Auxiliar/')
addpath('lib/poly2D/')
graf = 1;

%% Estimate displacement from consecutive sonograms

% Simulation paramters
% Siemens VF10-5 = [fs, c_c, f0, n_elem, elem_w, elv_focus, tx_focus] 
sim_data = [50e6, 1540, 6.67e6, 128, 0.3e-3, 10e-3, 10e-3]; 

% Import sonograms
E = 6000;    % Young's Modulus [Pa]
density = 5; % [scatters/mm^3]
load(sprintf('output/s_%d_d%d_3d.mat', E, density),...
    'sonograms', 'img_x', 'img_z');

if graf
    % Show first sonogram
    fig = figure(1);
    fig.Position = [0, 750, 750, 750];
    image_obj = imagesc(img_z - 1e-3, img_x, squeeze(sonograms(1, :, :)));
    ylabel('Lateral Distance (X) [m]')
    xlabel('Depth (Z) [m]')
    title('Sonogram with full FOV')
    axis equal;
    colormap(jet)
    colorbar
end

%% Estimate displacement
fine_prec = 0.05;
[est_x, est_z, u_x_est, u_z_est] = estimate_u(img_x, img_z, ...
                                               sonograms, fine_prec, graf);

%% Save results

save(sprintf('output/u_est_%d_d%d.mat', E, density),...
    'u_x_est', 'u_z_est', 'est_x', 'est_z')

%% Show results

% Import original displacement
E = 4000;    % Young's Modulus [Pa]
density = 5; % [scatters/mm^3]
 
% Load Real and estimated displacement
[t_dim, x_dim, y_dim, z_dim, u_x, u_z] = load_mid_u(sprintf("input/u_%d.h5", E));
load(sprintf('output/u_est_%d_d%d.mat', E, density),...
    'u_x_est', 'u_z_est', 'est_x', 'est_z')

if graf
  
    fig = figure(3);
    fig.Position = [500, 0, 500, 500];
    fig_mesh1 = mesh(z_dim, x_dim, squeeze(u_z(2, :, :)));
    title('Real Z Displacement')
    zlim([-1e-4, 2e-3])
    clim([-1e-4, 2e-3])
    view(120, 30)
    colormap("jet")
    colorbar
    
    filt_sz = 5;

    fig = figure(2);
    fig.Position = [0, 0, 500, 500];
    fig_mesh2 = mesh(est_z, est_x, conv2(squeeze(u_z_est(1, :, :)), ones(filt_sz)/filt_sz^2, 'same'));
    title('Estimated Z Displacement')
    zlim([-1e-4, 2e-3])
    clim([-1e-4, 2e-3])
    view(120, 30)
    colormap("jet")
    colorbar
    
    for t_i = 2:size(u_z_est, 1)
        pause(0.05)
        fig_mesh1.ZData = squeeze(u_z(t_i + 1, :, :));
        fig_mesh2.ZData = conv2(squeeze(u_z_est(t_i, :, :)), ones(filt_sz)/filt_sz^2, 'same');
    end
end
