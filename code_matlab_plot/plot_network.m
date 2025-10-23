
function [] = plot_network(x, neighbour_matrix, is_cross_linker, is_entangle, fig_name)
    cmap = colormap(jet);
    dof = size(x, 1);
    
    fig = figure; hold on; axis equal;

    for i = 1:dof
    for j = 1:4
        k = neighbour_matrix(i,j);
        if (i ~= k)
            line([x(i,1), x(k,1)], [x(i,2), x(k,2)], 'Color', cmap(1,:), 'linewidth', 2);
        end
    end
    end
    for i = 1:dof
        if (~is_entangle)
            plot(x(i,1), x(i,2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 4, 'linewidth', 1.2);
        else
            if (is_cross_linker(i) == 1)
                plot(x(i,1), x(i,2), 'o', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'MarkerSize', 4, 'linewidth', 1.2);
            else
                plot(x(i,1), x(i,2), 'x', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'MarkerSize', 4, 'LineWidth', 2);
            end
        end
    end
    
    % xlim([-2,82])
    
    % make the figure full screen
    set(fig, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);

    saveas(gcf, fig_name, 'fig');    
    saveas(gcf, fig_name, 'pdf');

    % close(fig);
end

