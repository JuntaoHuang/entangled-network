
function strech_info = compute_stretch_no_entangle(x, neighbour_matrix, x_init)

    dof = size(x, 1);
    strech_info = zeros(dof, 4);
    
    for i = 1:dof

        for j = 1:4
        
            k = neighbour_matrix(i,j);
            
            if (i ~= k)
                strech_info(i, j) = dist_point(x(i,:), x(k,:)) / dist_point(x_init(i,:), x_init(k,:));
            end
        end

    end    

end

