function [est_x, est_z, u_x_est, u_z_est] = estimate_u(img_x, img_z, ...
                                        sonograms, win_len, win_hop, ...
                                        fine_prec, graf)
%ESTIMATE_U Estimate X & Z displacements from a sequence of sonograms. Uses
%interpolation of the correlation function to achive sub-sample accuracy.

% Reduce safety gap
sonograms = sonograms(:, :, img_z >= 2e-3);
img_z = img_z(img_z >= 2e-3);

% Generate windows
dx = img_x(2) - img_x(1);
dz = img_z(2) - img_z(1);
w_x = 1 + 2 * ceil(win_len / dx / 2); % Uneven window length [samples]
w_z = 1 + 2 * ceil(win_len / dz / 2); % Uneven window length [samples]

% Prepare borders with reflection
sonograms = [flip(sonograms(:, 2:(w_x-1)/2, :), 2), sonograms];
img_x = [-flip(img_x(2:(w_x-1)/2)), img_x];

% Define window hop
hop_x = max(floor(win_hop / dx), 1); % Window hop size [samples]
hop_z = max(floor(win_hop / dz), 1); % Window hop size [samples]
N_x = floor((length(img_x) - w_x) / hop_x); % Number of windows
N_z = floor((length(img_z) - w_z) / hop_z); % Number of windows

% Create estimation and correlation dimensions
est_x = (0:N_x-1) * hop_x * dx;
est_z = ((0:N_z-1) * hop_z + ((w_z - 1) / 2)) * dz;

corr_x = (-w_x+1:w_x-1) * dx;
corr_z = (-w_z+1:w_z-1) * dz;

% Create coarse and fine grids
coarse_dim = -2:2;
[coarse_x, coarse_z] = meshgrid(coarse_dim, coarse_dim);
fine_dim = (-2:fine_prec:2);
[fine_x, fine_z] = meshgrid(fine_dim, fine_dim);

% Estimate displacement
u_x_est = zeros(size(sonograms, 1) - 1, N_x, N_z);
u_z_est = zeros(size(sonograms, 1) - 1, N_x, N_z);

%imagesc(img_z, img_x, squeeze(sonograms(1, :, :)));

for t = 1:size(sonograms, 1)-1
    if graf; tic(); end
    for x = 1:N_x
        for z = 1:N_z

            % Calculate correlation between windows
            win_x = (1:w_x) + (x - 1) * hop_x;
            win_z = (1:w_z) + (z - 1) * hop_z;

            win_corr = xcorr2(squeeze(sonograms(1, win_x, win_z)), ...
                              squeeze(sonograms(t + 1, win_x, win_z)));

            win_corr = win_corr / max(abs(win_corr), [], 'all');

            % Find maximum with coarse precission
            [~, max_corr] = max(win_corr, [], "all");
            [max_x, max_z] = ind2sub(size(win_corr), max_corr);

            % Establish maximum displacement limit
            max_x = min(max(max_x, 3), 2 * w_x - 3);
            max_z = min(max(max_z, 3), 2 * w_z - 3);

            % Fit correlation to polynomial around coarse maximum
            xcorr_p = polyFit2D(win_corr((-2:2) + max_x, (-2:2) + max_z), ...
                          coarse_x, coarse_z, 4, 4);
            
            % Evaluate polynomial with fine precission
            fine_xcorr = polyVal2D(xcorr_p, fine_x, fine_z, 4, 4);

            % Find maximum in fine grid
            [~, max_poly] = max(fine_xcorr, [], 'all');
            [max_dx, max_dz] = ind2sub(size(fine_xcorr), max_poly);

            % Calculate displacement
            u_x_est(t, x, z) = - corr_x(max_x) - fine_dim(max_dx) * dx;
            u_z_est(t, x, z) = - corr_z(max_z) - fine_dim(max_dz) * dz;

            if graf
            figure(1)
            subplot(1, 2, 1)
            imagesc(img_z(win_z), img_x(win_x), squeeze(sonograms(1, win_x, win_z)));
            title("Pre-deform Window")
            zlim([0, 1e-22])
            subplot(1, 2, 2)
            imagesc(img_z(win_z), img_x(win_x), squeeze(sonograms(t + 1, win_x, win_z)));
            title("Post-deform Window")
            zlim([0, 1e-22])

            figure(2)
            imagesc(corr_z, corr_x, win_corr);
            zlim([0, 1])
            title("Correlation between windows")
            
            figure(3)
            imagesc((-2:fine_prec:2), (-2:fine_prec:2), fine_xcorr);
            title("Fine-gridded correlation near coarse maximum")

            fprintf("Coarse : %.2g\nFine : %.2g", -corr_z(max_z), -fine_dim(max_dz) * dz)
            figure(1)
            waitforbuttonpress();
            end
        end
    end
    if graf
        it_time = toc();
        fprintf("\nProgreso = %.1f %%\nRestante = %.1f min\n", ...
                                 t / (size(sonograms, 1)-1) * 100, ...
                                 it_time * (size(sonograms, 1)-1 - t) / 60)
    end
end
end
