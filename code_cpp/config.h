#ifndef CONFIG_H
#define CONFIG_H

#include "utils.h"

void generate_entangle_network(
    int num_layer,
    int num_width,
    int max_chain_num,
    int max_node_num,
    std::vector<int>& chain_num,
    std::vector<std::vector<std::vector<int>>>& chain_info);

void generate_entangle_network_chain_info(
    int num_layer,
    int num_width,
    int max_chain_num,
    int max_node_num,
    const std::vector<int> & is_cross_linker,
    const std::vector<int> & entangle_orientation,
    const std::vector<std::vector<int>> & neighbour_matrix,
    int & real_max_node_num,
    std::vector<int>& chain_num,
    std::vector<std::vector<std::vector<int>>>& chain_info);

void compute_chain_init_len(
    int num_layer,
    int num_width,
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<int> & chain_num,
    const std::vector<std::vector<std::vector<int>>> & chain_info,
    std::vector<std::vector<double>> & dist_matrix_init,
    std::vector<std::vector<double>> & chain_init_len);

void add_defect(
    int num_layer,
    int num_width,
    const std::vector<int> & is_cross_linker,
    const std::vector<int> & entangle_orientation,
    std::vector<std::vector<int>>& neighbour_matrix
    );

void generate_network_neighbour_matrix(
    int num_layer,
    int num_width,
    const std::vector<std::vector<double>>& x,
    std::vector<std::vector<int>>& neighbour_matrix,
    int is_fracture);

void get_chain_node_list(const int current_node_i, const int current_node_j, 
        const int neighbour_node_i, const int neighbour_node_j,
        const int num_width,
        const std::vector<int> & entangle_orientation,
        const std::vector<int> & is_cross_linker,
        std::vector<int> & chain_node_list);

void get_neighbour_node_list(const int i, const int j, const std::vector<int> & neighbour_list,
    const int num_width,
    std::vector<int> & neighbour_node_i, std::vector<int> & neighbour_node_j);

void get_neighbour_node_orientation(const int current_node_i, const int current_node_j,
                                    const int neighbour_node_i, const int neighbour_node_j,
                                    std::string & neighbour_node_orientation);

void get_next_neighbour_node_coordinates(const int current_node_i, const int current_node_j,
                                         const int neighbour_node_i, const int neighbour_node_j,
                                         const int entangle_orientation_neighbour_node,
                                         int & next_neighbour_node_i, int & next_neighbour_node_j);

void initial_cross_linker_random(const int num_layer, const int num_width,
                                const double prob_cross_linker,
                                const int random_seed,
                                const int is_fracture,
                                std::vector<int> & is_cross_linker);

void initial_cross_linker_random_min_gap(const int num_layer, const int num_width,
                                const double prob_cross_linker,
                                const int random_seed,
                                const int min_gap_x,
                                const int min_gap_y,
                                const int is_fracture,
                                std::vector<int> & is_cross_linker);

void initial_cross_linker_random_min_gap_2(const int num_layer, const int num_width,
                                const double prob_cross_linker,
                                const int random_seed,
                                const int is_fracture,
                                std::vector<int> & is_cross_linker);

void initial_cross_linker_quasi_random(const int num_layer, const int num_width,
                                const double prob_cross_linker,
                                const int random_seed,
                                std::vector<int> & is_cross_linker);

double compute_cross_linker_ratio_inner_domain(int num_layer, int num_width, const std::vector<int> & is_cross_linker);

void initial_cross_linker_periodic_unit(const int num_layer, const int num_width,
                                        const int unit_size_x, const int unit_size_y,
                                        const int num_cross_linker_in_unit,
                                        const int is_fracture,
                                        std::vector<int> & is_cross_linker);

#endif