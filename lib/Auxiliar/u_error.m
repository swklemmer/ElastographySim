function [u_err, mean_err] = u_error(est_x, est_z, u_est, ...
                                            real_x, real_z, u_real)
%U_ERROR Calculates absolute error between estimated and real displacement.
% Both are 2D matrices.

% Interpolate real displacement in estimation grid
u_inter = interp2(real_z, real_x, u_real, est_z, est_x', 'linear');

% Calculate absolute error
u_err = abs(u_est - u_inter);

% Calculate mean error where displacement is significant (> 1%)
mean_err = mean(u_err(u_inter > 0.01), 'all');
end

