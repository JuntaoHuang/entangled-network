

function plot_network_stretch_no_entangle_stretch_max(x, neighbour_matrix, x_init, fig_name, stretch_max)
    
    strech_info = compute_stretch_no_entangle(x, neighbour_matrix, x_init);
    
    fig = figure; hold on; axis equal;
% set(fig, 'Visible', 'off');

    dof = size(x, 1);

ratio = 4;
    num_layer = round(sqrt(dof / ratio));
    num_width = num_layer * ratio;
    
    cmap = jet(256);
    
    % compute min and max stretch
stretch_min = 1.0;

    plotted_chain = [-1, -1];

    for i = 1:dof
        for l = 1:4
    
            j = neighbour_matrix(i,l);

            if (j~=i)

                node_pair_1 = [i, j];
                node_pair_2 = [j, i];

                if (~(any(ismember(plotted_chain, node_pair_1, 'rows')) ...
                        || any(ismember(plotted_chain, node_pair_2, 'rows'))))

                    plotted_chain(end+1, :) = node_pair_1;
                
                    data_x = [x(i, 1), x(j, 1)];
                    data_y = [x(i, 2), x(j, 2)];
                    data_c = strech_info(i,l);

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

    % make right boundary in black color
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
    c.Label.String = 'Value';

    % for i = 1:dof
    %     plot(x(i,1), x(i,2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 2.5, 'linewidth', 1.2);
    % end

% xlim([0, 80]);
% ylim([0, 100]);

    % make the figure full screen
    set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    
    saveas(gcf, fig_name, 'fig');
    saveas(gcf, fig_name, 'pdf');

    % close(fig);
end
