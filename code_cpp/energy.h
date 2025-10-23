#ifndef ENERGY_H
#define ENERGY_H

#include "utils.h"

double compute_total_energy_no_entangle(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<std::vector<double>>& dist_matrix_init);

double compute_total_energy_no_entangle_inner_domain(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const int num_layer, const int num_width,
    const std::vector<std::vector<double>>& dist_matrix_init);

double compute_total_energy_entangle(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<int>& chain_num,
    const std::vector<std::vector<double>>& chain_init_len,
    const std::vector<std::vector<std::vector<int>>>& chain_info);

double compute_total_energy_entangle_inner_domain(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<int>& chain_num,
    const std::vector<std::vector<double>>& chain_init_len,
    const std::vector<std::vector<std::vector<int>>>& chain_info,
    const int num_layer, const int num_width);

double compute_max_stretch_no_entangle(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<std::vector<double>>& dist_matrix_init);

double compute_max_stretch_no_entangle_inner_domain(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const int num_layer, const int num_width,
    const std::vector<std::vector<double>>& dist_matrix_init);

double compute_max_stretch_entangle(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<int>& chain_num,
    const std::vector<std::vector<double>>& chain_init_len,
    const std::vector<std::vector<std::vector<int>>>& chain_info);

double compute_max_stretch_entangle_inner_domain(
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<int>& chain_num,
    const std::vector<std::vector<double>>& chain_init_len,
    const std::vector<std::vector<std::vector<int>>>& chain_info,
    const int num_layer, const int num_width);

#endif