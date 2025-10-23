
function chain_init_len_list = compute_chain_init_len_distribution_entangle(chain_num, chain_init_len, chain_info, num_layer, num_width, is_cross_linker)

    [dof, ~, ~] = size(chain_info);

    chain_init_len_list = zeros(dof * 4, 1);
    
    chain_index = 1;

    for i_layer = 1:num_layer

        for i_width = 1:num_width
            
            i = i_width+(i_layer-1)*num_width;
            
            % only take the nodes in the center part into account
            if (i_width >= num_width*0.25 && i_width <= num_width*0.75 && i_layer >= num_layer*0.25 && i_layer <= num_layer*0.75 && is_cross_linker(i) == 1)

                for j = 1:chain_num(i)
                    
                    chain_init_len_list(chain_index) = chain_init_len(i,j);
                    chain_index = chain_index + 1;
                end
            end
        end
    end    

    
    chain_init_len_list(chain_index:end) = [];
end
