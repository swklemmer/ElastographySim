function [u_error, total_error] = u_error(est_x, est_z, u_est, ...
                                            real_x, real_z, u_real)
%U_ERROR Calculates absolute error between estimated and real displacement.
% Both are 2D matrices.

% Interpolate real displacement in estimation grid
u_inter = interp2(real_z, real_x, u_real, est_z, est_x', 'linear');

% Calculate absolute error
u_error = abs(u_est - u_inter);

% Calculate total error
total_error = sum(u_error, 'all');
end

