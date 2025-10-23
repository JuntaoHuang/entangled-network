
function [i, j] = global_to_local_index(node_index, num_width)

    j = ceil(node_index / num_width);
    i = node_index - (j - 1) * num_width;

end 