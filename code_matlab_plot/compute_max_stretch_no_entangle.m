

function max_stretch = compute_max_stretch_no_entangle(x, neighbour_matrix, x_init, num_layer, num_width)

    % dof = size(x, 1);
    max_stretch = 0.0;
    
    for j = 1:num_layer

        for i = 1:num_width
            
            node_index = i+(j-1)*num_width;

            % only take the nodes in the center part into account
            if (i >= num_width*0.25 && i <= num_width*0.75 && j >= num_layer*0.25 && j <= num_layer*0.75)
                
                for k = 1:4
                    neighbour_index = neighbour_matrix(node_index, k);
                    if (node_index ~= neighbour_index)
                        
                        stretch = dist_point(x(node_index,:), x(neighbour_index,:)) / dist_point(x_init(node_index,:), x_init(neighbour_index,:));
                        
                        if (stretch > max_stretch)
                            max_stretch = stretch;
                        end
                    end                
                end
                
            end

        end

    end

end

