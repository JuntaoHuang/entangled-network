
function d = get_dist_crack_tip(node_i, node_j, num_layer, num_width)

    crack_center_x = round(num_width / 2) + 0.5;
    crack_center_y = round(num_layer / 2) + 0.5;    

    [node_i_x, node_i_y] = global_to_local_index(node_i, num_width);
    [node_j_x, node_j_y] = global_to_local_index(node_j, num_width);

    chain_x = (node_i_x + node_j_x) / 2.0;
    chain_y = (node_i_y + node_j_y) / 2.0;

    dist_to_crack_center = dist_linf(chain_x, chain_y, crack_center_x, crack_center_y);
    
    d = floor(dist_to_crack_center) + 1;
end

% for two points p = (p_x, p_y) and p0 = (p0_x, p0_y)
% compute l_infty distance
% |p - p0|_inf = max(|p_x - p0_x|, |p_y - p0_y|)
function d = dist_linf(p_x, p_y, p0_x, p0_y)
   
    d = max(abs(p_x - p0_x), abs(p_y - p0_y));
    
end
