
function chain_index_i = get_chain_index_connect_two_node(node_i, node_j, chain_info, chain_num)
% this function is to find the chain index connecting two nodes
% 
% input:
% node_i: node index i
% node_j: node index j
% chain_info: chain information
% chain_num: number of chains for each node
% 
% output:
% chain_index_i: chain index for node_i

    chain_index_i = -1;

    for j = 1:chain_num(node_i)

        node_in_chain = reshape(chain_info(node_i, j ,:), 1, []);
        
        node_in_chain(node_in_chain==0) = [];

        if ismember(node_j, node_in_chain)
            chain_index_i = j;
            return;
        end
    end
    
    if chain_index_i < 0
        error('Error: cannot find the chain index connecting two nodes');
    end 
end

