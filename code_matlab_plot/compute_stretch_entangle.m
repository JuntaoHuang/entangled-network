
function strech_info = compute_stretch_entangle(x, neighbour_matrix, chain_num, chain_init_len, chain_info)

    dist_matrix = compute_dist_matrix(x, neighbour_matrix);

    [dof, max_chain_num, ~] = size(chain_info);

    strech_info = zeros(dof, max_chain_num);    
    
    for i = 1:dof
    
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
            strech_info(i,j) = chain_len/chain_init_len(i,j);            

        end    
    end

end

