function [fig_real, fig_est, fig_err] = show_disp_error(...
                                              real_x, real_z, u_real, ...
                                              est_x, est_z, u_est, ...
                                              u_err, mean_err, phan_dim)
%SHOW_DISPLACEMENT Show meshplot of real and estimated displacement

fig = tiledlayout(1, 3,'TileSpacing','Compact','Padding','Compact');
sgtitle('Normalized Z Displacement') 

ax1 = nexttile;
fig_real = mesh(real_z * 1e3, real_x * 1e3, u_real);
xlim([0 phan_dim(1) * 1e3 + 2])
ylim([0 phan_dim(2) * 1e3 + 2])
zlim([0, 1.2])
clim([0, 1.2])
view(30, 30)
title('Real')

ax2 = nexttile;
fig_est = mesh(est_z * 1e3, est_x * 1e3, u_est);
title('Estimated')
ylabel('Lateral Distance (X) [mm]')
xlabel('Depth (Z) [mm]')

ax3 = nexttile;
fig_err = mesh(est_z * 1e3, est_x * 1e3, u_err);
title(sprintf('Mean Abs. Error: %.2f A.U', mean_err))
colorbar

Link = linkprop([ax1, ax2, ax3],...
        {'CameraUpVector', 'CameraPosition', 'CameraTarget', ...
        'XLim', 'YLim', 'ZLim'});
setappdata(fig, 'StoreTheLink', Link);
end

