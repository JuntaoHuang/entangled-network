

function plot_network_energy_entangle_init(x, x_init, neighbour_matrix, chain_num, chain_init_len, chain_info, is_cross_linker, energy_function, fig_name)
    
    strech_info = compute_stretch_entangle(x, neighbour_matrix, chain_num, chain_init_len, chain_info);
    
    fig = figure; hold on; axis equal;
% set(fig, 'Visible', 'off');    
    dof = size(x, 1);
    
ratio = 4;
    num_layer = round(sqrt(dof / ratio));
    num_width = num_layer * ratio;

    cmap = jet(256);

    % compute min and max stretch
    stretch_min = 10.0;
    stretch_max = 0.0;
    for i = 1:dof    
        for j = 1:chain_num(i)
            stretch_min = min(stretch_min, strech_info(i,j));
            stretch_max = max(stretch_max, strech_info(i,j));
        end
    end
stretch_min = 1.0;
stretch_max = 5.0;
energy_min = energy_function(stretch_min);
energy_max = energy_function(stretch_max);

    plotted_chain = [-1, -1];

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
                    
                    data_x = [x_init(node_in_chain(k), 1); x_init(node_in_chain(k+1), 1)];
                    data_y = [x_init(node_in_chain(k), 2); x_init(node_in_chain(k+1), 2)];
                    data_c = strech_info(i,j);                    

if data_c <= stretch_min
    data_c = stretch_min;
end
if data_c >= stretch_max
    data_c = stretch_max;
end                
                    data_c = energy_function(data_c);

                    normalized_energy = (data_c - energy_min) / (energy_max - energy_min);
    
                    colorIndex = round(normalized_energy * (size(cmap, 1) - 1)) + 1;
                    plot(data_x, data_y, 'Color', cmap(colorIndex, :), 'LineWidth', 1);
                
                end
            end
        end
    end

    % make right boundary in black color
    index_i = num_width;
    for index_j = 1:(num_layer-1)
        index = (index_j - 1) * num_width + index_i;
        index_next = index_j * num_width + index_i;
        
        data_x = [x_init(index, 1); x_init(index_next, 1)];
        data_y = [x_init(index, 2); x_init(index_next, 2)];

        plot(data_x, data_y, 'k', 'LineWidth', 1);
    end
    
    % colormap('jet');    % assign the colormap
    % shading flat        % so each line segment has a plain color
    % view(2)             % set view in X-Y plane
    % colorbar
    % clim([4 5]);

    % Add a colorbar to show the mapping of color to original values
    colormap(cmap);
    c = colorbar;
    
    % Set colorbar ticks to reflect original values
    tickValues = linspace(0.0, 1.0, 6);
    c.Ticks = tickValues; % Normalize tick values
    c.TickLabels = arrayfun(@num2str, tickValues, 'UniformOutput', false); % Set original values as labels    
    % c.Label.String = 'Value';

    % for i = 1:dof
    %     if (is_cross_linker(i) == 1)
    %         plot(x_init(i,1), x_init(i,2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 2.5, 'linewidth', 1.2);
    %     else
    %         plot(x_init(i,1), x_init(i,2), 'x', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 2.5, 'LineWidth', 1.6);
    %     end
    % end

    % make the figure full screen
    set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    saveas(gcf, fig_name, 'fig');
    saveas(gcf, fig_name, 'pdf');

    % close(fig);
end
