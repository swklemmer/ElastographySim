function show_sonogram(img_x, img_z, sonograms, amp_lim, phan_dim, ...
                       time, title_txt)
%SHOW_SONOGRAM Plots image of sonogram.
imagesc(img_z * 1e3, img_x * 1e3, ...
        squeeze(sonograms(time, :, :)));
ylabel('Lateral Distance (X) [mm]')
xlabel('Depth (Z) [mm]')
title(title_txt)
axis equal
xlim([0 phan_dim(1) * 1e3 + 4])
ylim([0 phan_dim(2) * 1e3 + 2])
clim([0 amp_lim])
colormap(jet)
colorbar
end
