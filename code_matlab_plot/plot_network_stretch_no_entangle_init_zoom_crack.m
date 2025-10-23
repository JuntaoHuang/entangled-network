

function plot_network_stretch_no_entangle_init_zoom_crack(x, neighbour_matrix, x_init, fig_name)
    
    strech_info = compute_stretch_no_entangle(x, neighbour_matrix, x_init);
    
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
        for l = 1:4    
            j = neighbour_matrix(i,l);
            if (j~=i)
                stretch_min = min(stretch_min, strech_info(i,l));
                stretch_max = max(stretch_max, strech_info(i,l));
            end            
        end
    end
stretch_min = 1.0;
stretch_max = 5.0;

    %% plot chains near crack and stretch in color
    interval_in_x = 15;
    interval_in_y = 15;
    
    min_ix = round(num_width / 2) + 1 - interval_in_x;
    max_ix = round(num_width / 2) + 1 + interval_in_x;
    min_iy = round(num_layer / 2) - interval_in_y;
    max_iy = round(num_layer / 2) + 1 + interval_in_y;


    for ix = min_ix : max_ix
    for iy = min_iy : max_iy

        i = local_to_global_index(ix, iy, num_width);
        
        i = min(i, dof);
        i = max(i, 1);

        for l = 1:4
    
            j = neighbour_matrix(i,l);
            if (j~=i)

                node1 = i;
                node2 = j;

                [node1_ix, node1_iy] = global_to_local_index(node1, num_width);
                [node2_ix, node2_iy] = global_to_local_index(node2, num_width);

                if ((node1_ix <= max_ix) && (node1_ix >= min_ix) ...
                    && (node2_ix <= max_ix) && (node2_ix >= min_ix) ...
                    && (node1_iy <= max_iy) && (node1_iy >= min_iy) ...
                    && (node2_iy <= max_iy) && (node2_iy >= min_iy))
                
                    data_x = [x_init(i, 1); x_init(j, 1)];
                    data_y = [x_init(i, 2); x_init(j, 2)];
                    data_c = strech_info(i,l);

if data_c <= stretch_min
    data_c = stretch_min;
end
if data_c >= stretch_max
    data_c = stretch_max;
end

                    normalized_stretch = (data_c - stretch_min) / (stretch_max - stretch_min);
    
                    colorIndex = round(normalized_stretch * (size(cmap, 1) - 1)) + 1;
                    plot(data_x, data_y, 'Color', cmap(colorIndex, :), 'LineWidth', 8);
                
                end 
            end
        end
    end
    end

    % % make right boundary in black color
    % index_i = num_width;
    % for index_j = 1:(num_layer-1)
    %     index = (index_j - 1) * num_width + index_i;
    %     index_next = index_j * num_width + index_i;
    % 
    %     data_x = [x_init(index, 1); x_init(index_next, 1)];
    %     data_y = [x_init(index, 2); x_init(index_next, 2)];
    % 
    %     plot(data_x, data_y, 'k', 'LineWidth', 1);
    % end

    % Add a colorbar to show the mapping of color to original values
    colormap(cmap);
    c = colorbar;
    
    % Set colorbar ticks to reflect original values
    tickValues = linspace(1, stretch_max, 5); % Choose 5 evenly spaced ticks
    c.Ticks = (tickValues - stretch_min) / (stretch_max - stretch_min); % Normalize tick values
    c.TickLabels = arrayfun(@num2str, tickValues, 'UniformOutput', false); % Set original values as labels    
    % c.Label.String = 'Value'; 
    
    % for i = 1:dof
    %     plot(x_init(i,1), x_init(i,2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 2.5, 'linewidth', 1.2);        
    % end

    % make the figure full screen
    set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    
    saveas(gcf, fig_name, 'fig');
    saveas(gcf, fig_name, 'pdf');

    % close(fig);
end
