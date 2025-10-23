
function dist_matrix = compute_dist_matrix(x, neighbour_matrix)
    
    dof = size(x, 1);

    dist_matrix = sparse(dof, dof);
    for i = 1:dof
        for j = 1:4
            dist_matrix(i, neighbour_matrix(i,j)) = dist_point(x(i,:), x(neighbour_matrix(i,j),:));
        end
    end
    
end

