#include "energy.h"
#include "constitutive.h"

double compute_total_energy_no_entangle(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<std::vector<double>>& dist_matrix_init)
{
    int dof = x.size();
    double total_energy = 0.0;

    auto dist_matrix = compute_dist_matrix_simple(x, neighbour_matrix);

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            int neighbour_index = neighbour_matrix[i][j];
            if (i != neighbour_index)
            {
                double stretch = dist_matrix[i][j] / dist_matrix_init[i][j];
                double energy = energy_f(stretch) * dist_matrix_init[i][j];

                total_energy += energy / 2.0;
            }
        }
    }

    return total_energy;
}

double compute_total_energy_no_entangle_inner_domain(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const int num_layer, const int num_width,
    const std::vector<std::vector<double>>& dist_matrix_init)
{

    double total_energy = 0.0;

    auto dist_matrix = compute_dist_matrix_simple(x, neighbour_matrix);

    for (int j_layer = 0; j_layer < num_layer; ++j_layer)
    {
        for (int i_width = 0; i_width < num_width; ++i_width)
        {

            if (i_width >= 0.25 * num_width && i_width <= 0.75 * num_width && j_layer >= 0.25 * num_layer && j_layer <= 0.75 * num_layer)
            {
                int i = i_width + j_layer * num_width;

                for (int j = 0; j < 4; ++j)
                {
                    int neighbour_index = neighbour_matrix[i][j];
                    if (i != neighbour_index)
                    {
                        double stretch = dist_matrix[i][j] / dist_matrix_init[i][j];

                        double energy = energy_f(stretch) * dist_matrix_init[i][j];

                        total_energy += energy / 2.0;
                    }
                }
            }
        }
    }

    return total_energy;
}

double compute_total_energy_entangle(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<int>& chain_num,
    const std::vector<std::vector<int>>& chain_init_len,
    const std::vector<std::vector<std::vector<int>>>& chain_info)
{
    int dof = x.size();
    double total_energy = 0.0;

    auto dist_matrix = compute_dist_matrix_simple(x, neighbour_matrix);

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < chain_num[i]; ++j)
        {

            int num_node_in_chain = chain_info[i][j].size();

            double chain_len = 0.0;
            for (int k = 0; k < num_node_in_chain - 1; ++k)
            {
                int neighbor = chain_info[i][j][k + 1];
                for (int l = 0; l < 4; ++l)
                {
                    if (neighbor == neighbour_matrix[chain_info[i][j][k]][l])
                    {
                        chain_len += dist_matrix[chain_info[i][j][k]][l];
                        break;
                    }
                }
            }
            double stretch = chain_len / chain_init_len[i][j];

            double energy = energy_f(stretch) * chain_init_len[i][j];

            total_energy += energy / (int(chain_info[i][j].size()));
        }
    }

    return total_energy;    
}

double compute_total_energy_entangle_inner_domain(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<int>& chain_num,
    const std::vector<std::vector<double>>& chain_init_len,
    const std::vector<std::vector<std::vector<int>>>& chain_info,
    const int num_layer, const int num_width)
{

    double total_energy = 0.0;

    auto dist_matrix = compute_dist_matrix_simple(x, neighbour_matrix);

    for (int j_layer = 0; j_layer < num_layer; ++j_layer)
    {
        for (int i_width = 0; i_width < num_width; ++i_width)
        {

            if (i_width >= 0.25 * num_width && i_width <= 0.75 * num_width && j_layer >= 0.25 * num_layer && j_layer <= 0.75 * num_layer)
            {
                int i = i_width + j_layer * num_width;

                for (int j = 0; j < chain_num[i]; ++j)
                {

                    int num_node_in_chain = chain_info[i][j].size();

                    double chain_len = 0.0;
                    for (int k = 0; k < num_node_in_chain - 1; ++k)
                    {
                        int neighbor = chain_info[i][j][k + 1];
                        for (int l = 0; l < 4; ++l)
                        {
                            if (neighbor == neighbour_matrix[chain_info[i][j][k]][l])
                            {
                                chain_len += dist_matrix[chain_info[i][j][k]][l];
                                break;
                            }
                        }
                    }
                    double stretch = chain_len / chain_init_len[i][j];

                    double energy = energy_f(stretch) * chain_init_len[i][j];

                    total_energy += energy / (int(chain_info[i][j].size()));
                }
            }
        }
    }

    return total_energy;    
}

double compute_max_stretch_no_entangle(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<std::vector<double>>& dist_matrix_init)
{
    int dof = x.size();
    double max_stretch = 0.0;

    auto dist_matrix = compute_dist_matrix_simple(x, neighbour_matrix);

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            int neighbour_index = neighbour_matrix[i][j];
            if (i != neighbour_index)
            {
                double dist = dist_matrix[i][j];
                double stretch = dist / dist_matrix_init[i][j];

                max_stretch = std::max(max_stretch, stretch);
            }
        }
    }

    return max_stretch;
}

double compute_max_stretch_no_entangle_inner_domain(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const int num_layer, const int num_width,
    const std::vector<std::vector<double>>& dist_matrix_init)
{
    double max_stretch = 0.0;

    auto dist_matrix = compute_dist_matrix_simple(x, neighbour_matrix);

    for (int j_layer = 0; j_layer < num_layer; ++j_layer)
    {
        for (int i_width = 0; i_width < num_width; ++i_width)
        {

            if (i_width >= 0.25 * num_width && i_width <= 0.75 * num_width && j_layer >= 0.25 * num_layer && j_layer <= 0.75 * num_layer)
            {
                int i = i_width + j_layer * num_width;

                for (int j = 0; j < 4; ++j)
                {
                    int neighbour_index = neighbour_matrix[i][j];
                    if (i != neighbour_index)
                    {
                        double dist = dist_matrix[i][j];
                        double stretch = dist / dist_matrix_init[i][j];

                        max_stretch = std::max(max_stretch, stretch);
                    }
                }
            }
        }
    }

    return max_stretch;
}

double compute_max_stretch_entangle(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<int>& chain_num,
    const std::vector<std::vector<double>>& chain_init_len,
    const std::vector<std::vector<std::vector<int>>>& chain_info)
{
    int dof = x.size();
    double max_stretch = 0.0;

    auto dist_matrix = compute_dist_matrix_simple(x, neighbour_matrix);

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < chain_num[i]; ++j)
        {

            int num_node_in_chain = chain_info[i][j].size();

            double chain_len = 0.0;
            for (int k = 0; k < num_node_in_chain - 1; ++k)
            {
                int neighbor = chain_info[i][j][k + 1];
                for (int l = 0; l < 4; ++l)
                {
                    if (neighbor == neighbour_matrix[chain_info[i][j][k]][l])
                    {
                        chain_len += dist_matrix[chain_info[i][j][k]][l];
                        break;
                    }
                }
            }
            double stretch = chain_len / chain_init_len[i][j];

            max_stretch = std::max(max_stretch, stretch);
        }
    }

    return max_stretch;
}

double compute_max_stretch_entangle_inner_domain(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<int>& chain_num,
    const std::vector<std::vector<double>>& chain_init_len,
    const std::vector<std::vector<std::vector<int>>>& chain_info,
    const int num_layer, const int num_width)
{
    double max_stretch = 0.0;

    auto dist_matrix = compute_dist_matrix_simple(x, neighbour_matrix);

    for (int j_layer = 0; j_layer < num_layer; ++j_layer)
    {
        for (int i_width = 0; i_width < num_width; ++i_width)
        {

            if (i_width >= 0.25 * num_width && i_width <= 0.75 * num_width && j_layer >= 0.25 * num_layer && j_layer <= 0.75 * num_layer)
            {
                int i = i_width + j_layer * num_width;

                for (int j = 0; j < chain_num[i]; ++j)
                {

                    int num_node_in_chain = chain_info[i][j].size();

                    double chain_len = 0.0;
                    for (int k = 0; k < num_node_in_chain - 1; ++k)
                    {
                        int neighbor = chain_info[i][j][k + 1];
                        for (int l = 0; l < 4; ++l)
                        {
                            if (neighbor == neighbour_matrix[chain_info[i][j][k]][l])
                            {
                                chain_len += dist_matrix[chain_info[i][j][k]][l];
                                break;
                            }
                        }
                    }
                    double stretch = chain_len / chain_init_len[i][j];

                    max_stretch = std::max(max_stretch, stretch);
                }
            }

        }
    }

    return max_stretch;
}