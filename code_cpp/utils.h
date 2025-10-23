#ifndef UTILS_H
#define UTILS_H

#include <cmath>
#include <ctime>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip> 
#include <algorithm>
#include <limits>
#include <random>
#include <cassert>
#include <map>

int local_to_global_index(int i, int j, int num_width);

double dist_point(const std::vector<double>& p, const std::vector<double>& q);

double sum_abs_vec(const std::vector<double>& vec);

double sum_abs_mat(const std::vector<std::vector<double>>& mat);

std::vector<std::vector<double>> compute_dist_matrix(const std::vector<std::vector<double>>& x, const std::vector<std::vector<int>>& neighbour_matrix);

std::vector<std::vector<double>> compute_dist_matrix_simple(const std::vector<std::vector<double>>& x, const std::vector<std::vector<int>>& neighbour_matrix);

std::vector<int> remove_const(const std::vector<int>& vec, const int goal_const);

void save1DVectorDoubleToBinaryFile(
    const std::vector<double>& vec,
    const std::string& filename);

void save1DVectorIntToBinaryFile(
    const std::vector<int>& vec,
    const std::string& filename);

void save3DVectorDoubleToBinaryFile(
    const std::vector<std::vector<std::vector<double>>>& nestedVector,
    const std::string& filename);

void save3DVectorIntToBinaryFile(
    const std::vector<std::vector<std::vector<int>>>& nestedVector,
    const std::string& filename);

void save2DVectorIntToBinaryFile(
    const std::vector<std::vector<int>>& data,
    const std::string& filename);

void save2DVectorDoubleToBinaryFile(
    const std::vector<std::vector<double>>& matrix,
    const std::string& filename);

void removeNumber(std::vector<int>& vec, int numberToRemove);

double maxAbsoluteValue(const std::vector<std::vector<double>>& data);

double van_der_corput(const int n, const int base);

std::vector<std::vector<double>> generate_halton_seq(const int dim, const int n, const int skip);

std::vector<std::vector<int>> generate_halton_seq_2d_integer(const int num_points, const int x_lower_bound, const int x_upper_bound, const int y_lower_bound, const int y_upper_bound, const double rand_num);

void read_coordinate_from_file(const std::string& filename, std::vector<std::vector<double>>& x, int num_layer, int num_width);

#endif 