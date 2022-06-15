%% Importar
load("displacement_2.mat")

%% Animación

u = reshape(u, [size(u, 1), length(z_dim), length(x_dim)]);

fig = surf(x_dim, z_dim, squeeze(u(1, :, :)));
view(70, 45)
xlabel('Distancia axial (mm)')
ylabel('Distancia lateral (mm)')
zlabel('Desplazamiento (m)')
zlim([min(u(:,:,:), [], 'all'), max(u(:,:,:), [], 'all')]*1.1)

for t_i = 1:size(u, 1)
    fig.ZData = squeeze(u(t_i, :, :));
    pause(0.05)
end

