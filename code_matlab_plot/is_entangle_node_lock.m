
function is_lock = is_entangle_node_lock(x, entangle_node_i, is_cross_linker, chain_info, chain_num)

    is_lock = 0;

    if (is_cross_linker(entangle_node_i) == 1)
        return;
    end

    x_i = x(entangle_node_i, :);

    for j = 1:chain_num(entangle_node_i)

        node_in_chain = reshape(chain_info(entangle_node_i, j ,:), 1, []);
        
        node_in_chain(node_in_chain==0) = [];

        node_cross_link = [node_in_chain(1), node_in_chain(end)];

        for k = node_cross_link
            
            if (k ~= entangle_node_i)
                
                x_k = x(k, :);
                
                if (dist_point(x_i, x_k) <= 0.1)
                    is_lock = 1;
                end
            end
        end

    end

end

