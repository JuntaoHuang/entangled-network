#ifndef GRAIDENT_H
#define GRAIDENT_H

#include "utils.h"

void gradient_energy_no_entangle(std::vector<std::vector<double>>& grad,
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& xb,
    int num_width,
    int num_layer,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<std::vector<double>>& dist_matrix_init);

void gradient_energy_entangle(
    std::vector<std::vector<double>>& grad,
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<double>>& xb,
    const std::vector<std::vector<int>>& neighbour_matrix,
    int num_width,
    int num_layer,
    const std::vector<int>& chain_num,
    const std::vector<std::vector<double>>& chain_init_len,
    const std::vector<std::vector<std::vector<int>>>& chain_info);

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
    std::ofstream & file_output);

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
    const double stretch_crack);
    
#endif