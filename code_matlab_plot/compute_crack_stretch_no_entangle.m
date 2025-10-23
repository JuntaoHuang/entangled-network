
function crack_stretch = compute_crack_stretch_no_entangle(x, neighbour_matrix, x_init, num_layer, num_width)

    % lower node at crack tip
    i = round(num_width / 2) + 1;
    j = round(num_layer / 2);
    node_i = local_to_global_index(i, j, num_width);
    
    % upper node at crack tip
    i = round(num_width / 2) + 1;
    j = round(num_layer / 2) + 1;
    node_j = local_to_global_index(i, j, num_width);

    crack_stretch = compute_stretch_two_node_no_entangle(x, neighbour_matrix, x_init, node_i, node_j);

end