
function [strech_layer, num_chain_layer] = compute_chain_energy_layer_no_entangle(x, neighbour_matrix, x_init, num_layer, num_width)
    
    dof = size(x, 1);

    strech_info = compute_stretch_no_entangle(x, neighbour_matrix, x_init);

    strech_layer = zeros(round(num_width/2), 2*(num_layer + num_width));
    num_chain_layer = zeros(round(num_width/2), 1);
    
    for node_i = 1:dof
    
        for l = 1:4
        
            node_j = neighbour_matrix(node_i, l);
            
            if(node_i ~= node_j)
            
                d = get_dist_crack_tip(node_i, node_j, num_layer, num_width);

                num_chain_layer(d) = num_chain_layer(d) + 1;

                strech_layer(d, num_chain_layer(d)) = strech_info(node_i, l);
            
            end
        
        end
    end

end

