
function stretch = compute_stretch_two_node_no_entangle(x, neighbour_matrix, x_init, node_i, node_j)

    % get chain index connect two nodes
    chain_index_i = -1;
    for l = 1:4
        j = neighbour_matrix(node_i, l);
        if (j==node_j)
            chain_index_i = l;
        end
    end
                
    if (chain_index_i < 0)
        error('Error: cannot find the chain index connecting two nodes');
    end

    stretch = dist_point(x(node_i,:), x(node_j,:)) / dist_point(x_init(node_i,:), x_init(node_j,:));

end