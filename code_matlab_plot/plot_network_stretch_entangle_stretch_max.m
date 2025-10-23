

function plot_network_stretch_entangle_stretch_max(x, neighbour_matrix, chain_num, chain_init_len, chain_info, is_cross_linker, fig_name, stretch_max)
    
    strech_info = compute_stretch_entangle(x, neighbour_matrix, chain_num, chain_init_len, chain_info);
    
    fig = figure; hold on; axis equal;
% set(fig, 'Visible', 'off');
    dof = size(x, 1);
    
ratio = 4;
    num_layer = round(sqrt(dof / ratio));
    num_width = num_layer * ratio;

    cmap = jet(256);
    
    %% compute min and max stretch
stretch_min = 1.0;

    plotted_chain = [-1, -1];

    %% plot all the chains and stretch in color
    for i = 1:dof
    
        for j = 1:chain_num(i)
    
            % node index in this chain
            node_in_chain = reshape(chain_info(i,j,:), 1, []);
            
            % remove zero index
            node_in_chain(node_in_chain==0) = [];
    
            num_node_in_chain = length(node_in_chain);
            
            for k = 1:(num_node_in_chain-1)
                
                node_pair_1 = [node_in_chain(k), node_in_chain(k+1)];
                node_pair_2 = [node_in_chain(k+1), node_in_chain(k)];

                if (~(any(ismember(plotted_chain, node_pair_1, 'rows')) ...
                        || any(ismember(plotted_chain, node_pair_2, 'rows'))))

                    plotted_chain(end+1, :) = node_pair_1;
                
                    data_x = [x(node_in_chain(k), 1); x(node_in_chain(k+1), 1)];
                    data_y = [x(node_in_chain(k), 2); x(node_in_chain(k+1), 2)];
                    data_c = strech_info(i,j);

if data_c <= stretch_min
    data_c = stretch_min;
end
if data_c >= stretch_max
    data_c = stretch_max;
end
                    normalized_stretch = (data_c - stretch_min) / (stretch_max - stretch_min);    
                    colorIndex = round(normalized_stretch * (size(cmap, 1) - 1)) + 1;
                    plot(data_x, data_y, 'Color', cmap(colorIndex, :), 'LineWidth', 1);

                end
            end

        end
    end
    
    %% plot chain in the crack tip
    for crack_dist = 5:-1:0
        % lower node at crack tip
        i = round(num_width / 2) + 1;
        j = round(num_layer / 2);
        i = i + crack_dist;
        index_crack_lower = local_to_global_index(i, j, num_width);
    
        % upper node at crack tip
        i = round(num_width / 2) + 1;    
        j = round(num_layer / 2) + 1;
        i = i + crack_dist;
        index_crack_upper = local_to_global_index(i, j, num_width);
    
        chain_index_crack_lower = get_chain_index_connect_two_node(index_crack_lower, index_crack_upper, chain_info, chain_num);
    
        i = index_crack_lower;
        j = chain_index_crack_lower;
    
        node_in_chain = reshape(chain_info(i,j,:), 1, []);
    
        % remove zero index
        node_in_chain(node_in_chain==0) = [];
        num_node_in_chain = length(node_in_chain);
    
        for k = 1:(num_node_in_chain-1)
    
            data_x = [x(node_in_chain(k), 1); x(node_in_chain(k+1), 1)];
            data_y = [x(node_in_chain(k), 2); x(node_in_chain(k+1), 2)];
            data_c = strech_info(i,j);
    
            if data_c <= stretch_min
                data_c = stretch_min;
            end
            if data_c >= stretch_max
                data_c = stretch_max;
            end
    
            normalized_stretch = (data_c - stretch_min) / (stretch_max - stretch_min);
    
            colorIndex = round(normalized_stretch * (size(cmap, 1) - 1)) + 1;
            plot(data_x, data_y, 'Color', cmap(colorIndex, :), 'LineWidth', 1);
            
        end
    end

    % for k = 1:(num_node_in_chain)
    % 
    %     i = node_in_chain(k);
    % 
    %     if (is_cross_linker(i) == 1)
    %         plot(x(i,1), x(i,2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 2.5, 'linewidth', 1.2);
    %     else
    %         plot(x(i,1), x(i,2), 'x', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 2.5, 'LineWidth', 1.6);
    %     end
    % 
    % end
    
    %% make right boundary in black color
    index_i = num_width;
    for index_j = 1:(num_layer-1)
        index = (index_j - 1) * num_width + index_i;
        index_next = index_j * num_width + index_i;
        
        data_x = [x(index, 1); x(index_next, 1)];
        data_y = [x(index, 2); x(index_next, 2)];

        plot(data_x, data_y, 'k', 'LineWidth', 1);
    end

    % colormap('jet');    % assign the colormap
    % shading flat        % so each line segment has a plain color
    % view(2)             % set view in X-Y plane
    % colorbar

    % Add a colorbar to show the mapping of color to original values
    colormap(cmap);
    c = colorbar;

    % Set colorbar ticks to reflect original values
    tickValues = linspace(stretch_min, stretch_max, 5); % Choose 5 evenly spaced ticks
    c.Ticks = (tickValues - stretch_min) / (stretch_max - stretch_min); % Normalize tick values
    c.TickLabels = arrayfun(@num2str, tickValues, 'UniformOutput', false); % Set original values as labels    
    % c.Label.String = 'Value';

    for i = 1:dof
        if (is_cross_linker(i) == 1)
            % plot(x(i,1), x(i,2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 2.5, 'linewidth', 1.2);
            plot(x(i,1), x(i,2), 'b.');
        else
            % plot(x(i,1), x(i,2), 'x', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 2.5, 'LineWidth', 1.6);
            plot(x(i,1), x(i,2), 'r.');
        end
    end

% xlim([0, 80]);
% ylim([0, 100]);

    % make the figure full screen
    set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    
    saveas(gcf, fig_name, 'fig');
    saveas(gcf, fig_name, 'pdf');

    % close(fig);
end
