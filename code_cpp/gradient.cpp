#include "gradient.h"
#include "utils.h"
#include "constitutive.h"
#include "energy.h"

void gradient_energy_no_entangle(
    std::vector<std::vector<double>>& grad,
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& xb,
    int num_width,
    int num_layer,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<std::vector<double>>& dist_matrix_init)
{
    int dof = x.size();

    auto dist_matrix = compute_dist_matrix_simple(x, neighbour_matrix);

    grad.assign(dof, std::vector<double>(2, 0.0));

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            int neighbour_index = neighbour_matrix[i][j];
            if (i != neighbour_index)
            {
                double dist = dist_matrix[i][j];
                double stretch = dist / dist_matrix_init[i][j];

                grad[i][0] += f(stretch) / dist * (x[i][0] - x[neighbour_index][0]);
                grad[i][1] += f(stretch) / dist * (x[i][1] - x[neighbour_index][1]);
            }
        }
    }

    for (int i = 0; i < num_width; ++i)
    {
        for (int d = 0; d < 2; ++d)
        {
            grad[i][d] = x[i][d] - xb[i][d];
        }
    }

    for (int i = 0; i < num_width; ++i)
    {
        for (int d = 0; d < 2; ++d)
        {
            grad[dof - i - 1][d] = x[dof - i - 1][d] - xb[dof - i - 1][d];
        }
    }

    for (int i = 0; i < num_layer; ++i)
    {
        grad[num_width * i][0] = x[num_width * i][0] - xb[num_width * i][0];
    }

    for (int i = 0; i < num_layer; ++i)
    {
        grad[num_width * (i + 1) - 1][0] = x[num_width * (i + 1) - 1][0] - xb[num_width * (i + 1) - 1][0];
    }
}

void gradient_energy_entangle(
    std::vector<std::vector<double>>& grad,
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& xb,
    const std::vector<std::vector<int>>& neighbour_matrix,
    int num_width,
    int num_layer,
    const std::vector<int>& chain_num,
    const std::vector<std::vector<double>>& chain_init_len,
    const std::vector<std::vector<std::vector<int>>>& chain_info)
{
    int dof = x.size();

    auto dist_matrix = compute_dist_matrix_simple(x, neighbour_matrix);

    grad.assign(dof, std::vector<double>(2, 0.0));

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

            for (int k = 0; k < num_node_in_chain; ++k)
            {
                int node_in_chain = chain_info[i][j][k];
                for (int l = 0; l < 4; ++l)
                {
                    if (node_in_chain == neighbour_matrix[i][l] && node_in_chain != i)
                    {
                        grad[i][0] += f(stretch) * (x[i][0] - x[node_in_chain][0]) / dist_matrix[i][l];
                        grad[i][1] += f(stretch) * (x[i][1] - x[node_in_chain][1]) / dist_matrix[i][l];
                    }
                }
            }
        }
    }

    for (int i = 0; i < num_width; ++i)    
    {
        for (int d = 0; d < 2; ++d)
        {
            grad[i][d] = x[i][d] - xb[i][d];
        }
    }

    for (int i = 0; i < num_width; ++i)
    {
        for (int d = 0; d < 2; ++d)
        {
            grad[dof - i - 1][d] = x[dof - i - 1][d] - xb[dof - i - 1][d];
        }
    }

    for (int i = 0; i < num_layer; ++i)
    {
        grad[num_width * i][0] = x[num_width * i][0] - xb[num_width * i][0];
    }

    for (int i = 0; i < num_layer; ++i)
    {
        grad[num_width * (i + 1) - 1][0] = x[num_width * (i + 1) - 1][0] - xb[num_width * (i + 1) - 1][0];
    }
}

std::vector<std::vector<double>> solve_optimization_elasticity(const std::vector<std::vector<double>> & x_init,
    const std::vector<std::vector<double>> & xb,
    const bool is_entangle,
    const double dt_init,
    const int num_decay_dt,    
    const int max_time_step,
    const int print_time_interval,
    const int num_width, const int num_layer,
    const double error_tol,
    const std::vector<std::vector<int>> & neighbour_matrix,
    const std::vector<std::vector<double>> & dist_matrix_init,
    const std::vector<int> & chain_num,
    const std::vector<std::vector<double>> & chain_init_len,
    const std::vector<std::vector<std::vector<int>>> & chain_info,
    std::ofstream & file_output)
{   
    if (is_entangle)
    {
        file_output << "solve network with entangle \n";
    }
    else
    {
        file_output << "solve network without entangle \n";
    }

    const int dof = x_init.size();
    std::vector<std::vector<double>> x = x_init;

    bool is_solved = false;
    for (int i_decay_dt = 0; i_decay_dt < num_decay_dt; ++i_decay_dt)
    {
        double dt = dt_init * std::pow(0.1, i_decay_dt);

        file_output << "solve using time step " << dt << std::endl;
        file_output << "---------- \n";

        std::vector<std::vector<double>> grad(dof, std::vector<double>(2, 0.0));    
        std::vector<std::vector<double>> vel(dof, std::vector<double>(2, 0.0));     
        double beta = 0.9;  

        std::vector<double> total_energy_list_in_iteration(max_time_step/print_time_interval, 0.0);

        std::vector<std::vector<double>> x_previous_iteration(dof, std::vector<double>(2, 0.0));
        std::vector<std::vector<double>> x_current_iteration(dof, std::vector<double>(2, 0.0));
        for (int num_time_step = 0; num_time_step < max_time_step; ++num_time_step)
        {
            std::vector<std::vector<double>> grad(dof, std::vector<double>(2, 0.0));
            if (is_entangle)
            {
                gradient_energy_entangle(grad, x, xb, neighbour_matrix, num_width, num_layer, chain_num, chain_init_len, chain_info);
            }
            else
            {
                gradient_energy_no_entangle(grad, x, xb, num_width, num_layer, neighbour_matrix, dist_matrix_init);
            }

            for (int i = 0; i < dof; ++i)
            {
                vel[i][0] = beta * vel[i][0] + (1 - beta) * grad[i][0];
                vel[i][1] = beta * vel[i][1] + (1 - beta) * grad[i][1];
            }

            for (int i = 0; i < dof; ++i)
            {
                x[i][0] -= dt * vel[i][0];
                x[i][1] -= dt * vel[i][1];
            }

            if (num_time_step % print_time_interval == 0)
            {
                int time_index = num_time_step / print_time_interval;

                double l1_norm = sum_abs_mat(grad) / (grad.size() * grad[0].size());
                double linf_norm = maxAbsoluteValue(grad);

                double total_energy = 0.0;
                double max_stretch = 0.0;
                if (is_entangle)
                {
                    total_energy = compute_total_energy_entangle_inner_domain(x, neighbour_matrix, chain_num, chain_init_len, chain_info, num_layer, num_width);
                    max_stretch = compute_max_stretch_entangle_inner_domain(x, neighbour_matrix, chain_num, chain_init_len, chain_info, num_layer, num_width);
                }
                else
                {
                    total_energy = compute_total_energy_no_entangle_inner_domain(x, neighbour_matrix, num_layer, num_width, dist_matrix_init);
                    max_stretch = compute_max_stretch_no_entangle_inner_domain(x, neighbour_matrix, num_layer, num_width, dist_matrix_init);
                }

                total_energy_list_in_iteration[time_index] = total_energy;

                x_previous_iteration = x_current_iteration;
                x_current_iteration = x;

                file_output << std::scientific << std::setprecision(2)
                    << "num of time step " << num_time_step
                    << "; l1 norm of grad " << l1_norm
                    << "; linf norm of grad " << linf_norm
                    << "; total energy " << total_energy
                    << "; max stretch " << max_stretch
                    << std::endl;

                if (std::isnan(l1_norm) || std::isinf(l1_norm) || std::isnan(linf_norm) || std::isinf(linf_norm))
                {
                    file_output << "network generate nan or inf; use smaller time step for iteration. \n";
                    file_output << "---------- \n";
                    x = x_init;
                    break;
                }

                if (linf_norm < error_tol)
                {
                    is_solved = true;
                    file_output << "error tolerance for gradient is satisfied; network is solved. \n";
                    file_output << "---------- \n";
                    break;
                }

                if (time_index >= 1)
                {
                    if (total_energy_list_in_iteration[time_index] > total_energy_list_in_iteration[time_index - 1])
                    {
                        if (i_decay_dt == num_decay_dt - 1)
                        {

                            x = x_previous_iteration;
                            file_output << "energy is increasing; stop iteration. \n";
                            file_output << "---------- \n";
                            break;
                        }
                        else
                        {

                            x = x_previous_iteration;
                            file_output << "energy is increasing; use smaller time step for iteration. \n";
                            file_output << "---------- \n";
                            break;
                        }
                    }
                }    
            }

            if ((num_time_step == max_time_step - 1) && (i_decay_dt == num_decay_dt - 1))
            {
                double l1_norm = sum_abs_mat(grad) / (grad.size() * grad[0].size());
                double linf_norm = maxAbsoluteValue(grad);

                if (!(std::isnan(l1_norm) || std::isinf(l1_norm) || std::isnan(linf_norm) || std::isinf(linf_norm)))
                {
                    is_solved = true;
                    file_output << "iteration results (in the last and smallest time step) are not nan or inf; network is solved. \n";
                    file_output << "---------- \n";
                }
            }
        }

        if (is_solved)
        {
            break;
        }
    }

    return x;
}

std::vector<std::vector<double>> solve_optimization_fracture(const std::vector<std::vector<double>> & x_init,
    const std::vector<std::vector<double>> & xb,
    const bool is_entangle,
    const double dt_init,
    const int num_decay_dt,    
    const int max_time_step,
    const int print_time_interval,
    const int num_width, const int num_layer,
    const double error_tol,
    const std::vector<std::vector<int>> & neighbour_matrix,
    const std::vector<std::vector<double>> & dist_matrix_init,
    const std::vector<int> & chain_num,
    const std::vector<std::vector<double>> & chain_init_len,
    const std::vector<std::vector<std::vector<int>>> & chain_info,
    std::ofstream & file_output,
    const int is_stop_iter_max_stretch,
    const double stretch_crack)
{   
    if (is_entangle)
    {
        file_output << "solve network with entangle \n";
    }
    else
    {
        file_output << "solve network without entangle \n";
    }

    const int dof = x_init.size();
    std::vector<std::vector<double>> x = x_init;

    bool is_solved = false;
    for (int i_decay_dt = 0; i_decay_dt < num_decay_dt; ++i_decay_dt)
    {
        double dt = dt_init * std::pow(0.1, i_decay_dt);

        file_output << "solve using time step " << dt << std::endl;
        file_output << "---------- \n";

        std::vector<std::vector<double>> grad(dof, std::vector<double>(2, 0.0));    
        std::vector<std::vector<double>> vel(dof, std::vector<double>(2, 0.0));     
        double beta = 0.9;  

        std::vector<double> total_energy_list_in_iteration(max_time_step/print_time_interval, 0.0);
        std::vector<double> l1_norm_list_in_iteration(max_time_step/print_time_interval, 0.0);
        std::vector<double> linf_norm_list_in_iteration(max_time_step/print_time_interval, 0.0);

        std::vector<std::vector<double>> x_previous_iteration(dof, std::vector<double>(2, 0.0));
        std::vector<std::vector<double>> x_current_iteration(dof, std::vector<double>(2, 0.0));
        for (int num_time_step = 0; num_time_step < max_time_step; ++num_time_step)
        {
            std::vector<std::vector<double>> grad(dof, std::vector<double>(2, 0.0));
            if (is_entangle)
            {
                gradient_energy_entangle(grad, x, xb, neighbour_matrix, num_width, num_layer, chain_num, chain_init_len, chain_info);
            }
            else
            {
                gradient_energy_no_entangle(grad, x, xb, num_width, num_layer, neighbour_matrix, dist_matrix_init);
            }

            for (int i = 0; i < dof; ++i)
            {
                vel[i][0] = beta * vel[i][0] + (1 - beta) * grad[i][0];
                vel[i][1] = beta * vel[i][1] + (1 - beta) * grad[i][1];
            }

            for (int i = 0; i < dof; ++i)
            {
                x[i][0] -= dt * vel[i][0];
                x[i][1] -= dt * vel[i][1];
            }

            if (num_time_step % print_time_interval == 0)
            {
                int time_index = num_time_step / print_time_interval;

                double l1_norm = sum_abs_mat(grad) / (grad.size() * grad[0].size());
                double linf_norm = maxAbsoluteValue(grad);

                double total_energy = 0.0;
                double max_stretch = 0.0;
                if (is_entangle)
                {
                    total_energy = compute_total_energy_entangle_inner_domain(x, neighbour_matrix, chain_num, chain_init_len, chain_info, num_layer, num_width);
                    max_stretch = compute_max_stretch_entangle_inner_domain(x, neighbour_matrix, chain_num, chain_init_len, chain_info, num_layer, num_width);
                }
                else
                {
                    total_energy = compute_total_energy_no_entangle_inner_domain(x, neighbour_matrix, num_layer, num_width, dist_matrix_init);
                    max_stretch = compute_max_stretch_no_entangle_inner_domain(x, neighbour_matrix, num_layer, num_width, dist_matrix_init);
                }

                total_energy_list_in_iteration[time_index] = total_energy;
                l1_norm_list_in_iteration[time_index] = l1_norm;
                linf_norm_list_in_iteration[time_index] = linf_norm;

                x_previous_iteration = x_current_iteration;
                x_current_iteration = x;

                file_output << std::scientific << std::setprecision(2)
                    << "num of time step " << num_time_step
                    << "; l1 norm of grad " << l1_norm
                    << "; linf norm of grad " << linf_norm
                    << "; total energy " << total_energy
                    << "; max stretch " << max_stretch
                    << std::endl;

                if (std::isnan(l1_norm) || std::isinf(l1_norm) || std::isnan(linf_norm) || std::isinf(linf_norm))
                {
                    file_output << "network generate nan or inf; use smaller time step for iteration. \n";
                    file_output << "---------- \n";
                    x = x_init;
                    break;
                }

                if (linf_norm < error_tol)
                {
                    is_solved = true;
                    file_output << "error tolerance for gradient is satisfied; network is solved. \n";
                    file_output << "---------- \n";
                    return x;
                }

                if (time_index >= 1)
                {
                    if ((l1_norm_list_in_iteration[time_index] > l1_norm_list_in_iteration[time_index - 1]) && (total_energy_list_in_iteration[time_index] > total_energy_list_in_iteration[time_index - 1]))
                    {
                        if (i_decay_dt == num_decay_dt - 1)
                        {

                            x = x_previous_iteration;
                            file_output << "energy is increasing; stop iteration. \n";
                            file_output << "---------- \n";
                            break;
                        }
                        else
                        {

                            x = x_previous_iteration;
                            file_output << "energy is increasing; use smaller time step for iteration. \n";
                            file_output << "---------- \n";
                            break;
                        }
                    }
                }    

                if (is_stop_iter_max_stretch == 1)
                {                
                    if ((time_index >= 10) && (max_stretch > stretch_crack + 1.0))
                    {
                        file_output << "in fracture, max stretch " << max_stretch << " is greater than the stretch for crack " << stretch_crack << "; decrease the stretch for crack. \n";
                        file_output << "---------- \n";
                        return x;
                    }
                }
            }

            if ((num_time_step == max_time_step - 1) && (i_decay_dt == num_decay_dt - 1))
            {
                double l1_norm = sum_abs_mat(grad) / (grad.size() * grad[0].size());
                double linf_norm = maxAbsoluteValue(grad);

                if (!(std::isnan(l1_norm) || std::isinf(l1_norm) || std::isnan(linf_norm) || std::isinf(linf_norm)))
                {
                    is_solved = true;
                    file_output << "iteration results (in the last and smallest time step) are not nan or inf; network is solved. \n";
                    file_output << "---------- \n";
                    return x;
                }
            }
        }

        if (is_solved)
        {
            break;
        }
    }

    return x;
}