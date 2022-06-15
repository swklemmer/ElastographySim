function img_x = lateral_positions(sim_data, x_dim)
%LATERAL_POSITIONS Retrieve x positions from elements that are in front of
%phantom. Assumes x symmetry.

n_elem = sim_data(4);
elem_w = sim_data(5);

elem_x = ((1:n_elem) - (n_elem + 1) / 2) * elem_w;

% Find first and last element
x_0 = find(elem_x <= 0, 1, 'last');
x_f = find(elem_x > x_dim(end), 1, 'first');

img_x = elem_x(x_0:x_f);

end

