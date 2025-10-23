
function [crack_chain_len, crack_stretch] = compute_stretch_chain_len_two_node_entangle(x, chain_num, chain_init_len, chain_info, node_i, node_j)

    % strech_info = compute_stretch_entangle(x, neighbour_matrix, chain_num, chain_init_len, chain_info);

    % dof = size(x, 1);
    
    % num_layer = round(sqrt(dof / ratio_width_layer));
    % num_width = num_layer * ratio_width_layer;

    % % lower node at crack tip
    % i = round(num_width / 2) + 1;
    % j = round(num_layer / 2);
    % node_i = local_to_global_index(i, j, num_width);
    % 
    % % upper node at crack tip
    % i = round(num_width / 2) + 1;
    % j = round(num_layer / 2) + 1;
    % node_j = local_to_global_index(i, j, num_width);

    % get chain index connect two nodes at crack tip
    chain_index_i = get_chain_index_connect_two_node(node_i, node_j, chain_info, chain_num);

    % node index in this chain
    node_in_chain = reshape(chain_info(node_i, chain_index_i, :), 1, []);
    node_in_chain(node_in_chain==0) = [];

    num_node_in_chain = length(node_in_chain);
    crack_chain_len = num_node_in_chain - 1;

    % crack_stretch = strech_info(node_i, chain_index_i);

    chain_len_after_deform = 0;
    for k = 1:(num_node_in_chain-1)
        chain_len_after_deform = chain_len_after_deform + dist_point(x(node_in_chain(k),:), x(node_in_chain(k+1),:));
    end

    crack_stretch = chain_len_after_deform/chain_init_len(node_i, chain_index_i);
    
end