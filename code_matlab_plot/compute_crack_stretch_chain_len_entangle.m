
function [crack_chain_len, crack_stretch] = compute_crack_stretch_chain_len_entangle(x, chain_num, chain_init_len, chain_info, num_layer, num_width)

    % lower node at crack tip
    i = round(num_width / 2) + 1;
    j = round(num_layer / 2);
    node_i = local_to_global_index(i, j, num_width);
    
    % upper node at crack tip
    i = round(num_width / 2) + 1;
    j = round(num_layer / 2) + 1;
    node_j = local_to_global_index(i, j, num_width);

    [crack_chain_len, crack_stretch] = compute_stretch_chain_len_two_node_entangle(x, chain_num, chain_init_len, chain_info, node_i, node_j);
    
end