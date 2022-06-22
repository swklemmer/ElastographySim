function [u_err, mean_err] = estimation_error(real_x, real_z, u_real,...
                                            est_x, est_z, u_est, phan_dim)
%U_ERROR Calculates absolute error between estimated and real displacement.
% Both are 2D matrices.

% Interpolate real displacement in estimation grid
u_inter = interp2(real_z, real_x, u_real, est_z, est_x', 'linear');

% Calculate absolute error
u_err = abs(u_est - u_inter);

% Calculate mean error inside phantom

mean_err = mean(u_err(est_x >= 0 & est_x <= phan_dim(1),...
                      est_z >= 0 & est_z <= phan_dim(2)), 'all');
end

