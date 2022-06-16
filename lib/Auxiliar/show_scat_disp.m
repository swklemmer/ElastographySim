function show_scat_disp(scat_pos, new_scat_pos, phant_dim)
%SHOW_SCATTERERS Shows scatterplot of scatterers before and after
% displacement
scatter(scat_pos(:, 3) * 1e3, scat_pos(:, 1) * 1e3)
hold on
scatter(new_scat_pos(:, 3) * 1e3, new_scat_pos(:, 1) * 1e3)
hold off
axis ij
grid on
title('Scatterer Position')
axis equal
xlim([0 phant_dim(1) * 1e3 + 4])
ylim([0 phant_dim(2) * 1e3 + 2])
xlabel('Depth (Z) [mm]')
ylabel('Lateral Distance (X) [mm]')
legend({'Pre', 'Post'})
end

