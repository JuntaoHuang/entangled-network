
function chain_energy_list = compute_chain_energy_distribution_entangle(x, neighbour_matrix, chain_num, chain_init_len, chain_info, num_layer, num_width, is_cross_linker)

    dist_matrix = compute_dist_matrix(x, neighbour_matrix);

    [dof, ~, ~] = size(chain_info);

    chain_energy_list = zeros(dof * 4, 1);
    
    chain_index = 1;

    for i_layer = 1:num_layer

        for i_width = 1:num_width
            
            i = i_width+(i_layer-1)*num_width;
            
            % only take the nodes in the center part into account
            if (i_width >= num_width*0.25 && i_width <= num_width*0.75 && i_layer >= num_layer*0.25 && i_layer <= num_layer*0.75 && is_cross_linker(i) == 1)

                for j = 1:chain_num(i)

                    % node index in this chain
                    node_in_chain = reshape(chain_info(i,j,:), 1, []);
                    
                    % remove zero index
                    node_in_chain(node_in_chain==0) = [];
        
                    num_node_in_chain = length(node_in_chain);
        
                    % compute chain length after deformation
                    chain_len = 0;
                    for k = 1:(num_node_in_chain-1)
                        chain_len = chain_len + dist_matrix(node_in_chain(k), node_in_chain(k+1));
                    end
                    stretch = chain_len/chain_init_len(i,j);
                    
                    energy = energy_func(stretch) * chain_init_len(i,j);
                    
                    chain_energy_list(chain_index) = energy;
                    chain_index = chain_index + 1;
                end
            end
        end
    end    

    
    chain_energy_list(chain_index:end) = [];
end
