#include "utils.h"

int local_to_global_index(int i, int j, int num_width)
{
    return i + j * num_width;
}

double dist_point(const std::vector<double>& p, const std::vector<double>& q)
{
    return std::sqrt(std::pow(p[0] - q[0], 2) + std::pow(p[1] - q[1], 2));
}

double sum_abs_vec(const std::vector<double>& vec)
{
    double sum = 0.0;
    for (double val : vec)
    {
        if (std::isnan(val) || std::isinf(val))
        {
            sum = std::numeric_limits<double>::infinity();
            return sum;
        }
        sum += std::abs(val);
    }
    return sum;
}

double sum_abs_mat(const std::vector<std::vector<double>>& mat)
{
    double sum = 0.0;
    for (const auto& vec : mat)
    {
        sum += sum_abs_vec(vec);
    }
    return sum;
}

std::vector<std::vector<double>> compute_dist_matrix(const std::vector<std::vector<double>>& x, const std::vector<std::vector<int>>& neighbour_matrix)
{
    int dof = x.size();
    std::vector<std::vector<double>> dist_matrix(dof, std::vector<double>(dof, 0.0));
    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            int neighbour_index = neighbour_matrix[i][j];
            dist_matrix[i][neighbour_index] = dist_point(x[i], x[neighbour_index]);
        }
    }
    return dist_matrix;
}

std::vector<std::vector<double>> compute_dist_matrix_simple(const std::vector<std::vector<double>>& x, const std::vector<std::vector<int>>& neighbour_matrix)
{
    int dof = x.size();

    std::vector<std::vector<double>> dist_matrix(dof, std::vector<double>(4, 0.0));

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            int neighbour_index = neighbour_matrix[i][j];
            dist_matrix[i][j] = dist_point(x[i], x[neighbour_index]);
        }
    }

    return dist_matrix;
}

std::vector<int> remove_const(const std::vector<int>& vec, const int goal_const)
{
    std::vector<int> result;
    for (int val : vec)
    {
        if (val != goal_const)
        {
            result.push_back(val);
        }
    }
    return result;
}

void save3DVectorDoubleToBinaryFile(
    const std::vector<std::vector<std::vector<double>>>& data,
    const std::string& filename)
{
    std::ofstream outFile(filename, std::ios::binary);

    if (!outFile) {
        std::cerr << "Could not open the file for writing." << std::endl;
        return;
    }

    int x_dim = data.size();
    int y_dim = data[0].size();
    int z_dim = data[0][0].size();

    outFile.write(reinterpret_cast<const char*>(&x_dim), sizeof(int));
    outFile.write(reinterpret_cast<const char*>(&y_dim), sizeof(int));
    outFile.write(reinterpret_cast<const char*>(&z_dim), sizeof(int));

    for (const auto& mat : data) {
        for (const auto& row : mat) {
            outFile.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(double));
        }
    }

    outFile.close();
}

void save3DVectorIntToBinaryFile(
    const std::vector<std::vector<std::vector<int>>>& data,
    const std::string& filename)
{
    std::ofstream outFile(filename, std::ios::binary);

    if (!outFile) {
        std::cerr << "Could not open the file for writing." << std::endl;
        return;
    }

    int x_dim = data.size();
    int y_dim = data[0].size();
    int z_dim = data[0][0].size();

    outFile.write(reinterpret_cast<const char*>(&x_dim), sizeof(int));
    outFile.write(reinterpret_cast<const char*>(&y_dim), sizeof(int));
    outFile.write(reinterpret_cast<const char*>(&z_dim), sizeof(int));

    for (const auto& mat : data) {
        for (const auto& row : mat) {
            outFile.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(int));
        }
    }

    outFile.close();
}

void removeNumber(std::vector<int>& vec, int numberToRemove)
{
    auto newEnd = std::remove(vec.begin(), vec.end(), numberToRemove);

    vec.erase(newEnd, vec.end());
}

double maxAbsoluteValue(const std::vector<std::vector<double>>& data)
{

    double max_abs_value = 0.0;

    for (const auto& row : data)
    {
        for (const auto& value : row)
        {
            if (std::isnan(value) || std::isinf(value))
            {
                max_abs_value = std::numeric_limits<double>::infinity();
                return max_abs_value;
            }
            max_abs_value = std::max(max_abs_value, std::abs(value));
        }
    }

    return max_abs_value;
}

void save2DVectorIntToBinaryFile(
    const std::vector<std::vector<int>>& data,
    const std::string& filename)
{
    std::ofstream outFile(filename, std::ios::binary);

    if (!outFile) {
        std::cerr << "Could not open the file for writing." << std::endl;
        return;
    }

    int rows = data.size();
    int cols = data[0].size();

    outFile.write(reinterpret_cast<const char*>(&rows), sizeof(int));
    outFile.write(reinterpret_cast<const char*>(&cols), sizeof(int));

    for (const auto& row : data) {
        outFile.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(int));
    }

    outFile.close();
}

void save2DVectorDoubleToBinaryFile(
    const std::vector<std::vector<double>>& data,
    const std::string& filename)
{
    std::ofstream outFile(filename, std::ios::binary);

    if (!outFile) {
        std::cerr << "Could not open the file for writing." << std::endl;
        return;
    }

    int rows = data.size();
    int cols = data[0].size();

    outFile.write(reinterpret_cast<const char*>(&rows), sizeof(int));
    outFile.write(reinterpret_cast<const char*>(&cols), sizeof(int));

    for (const auto& row : data) {
        outFile.write(reinterpret_cast<const char*>(row.data()), row.size() * sizeof(double));
    }

    outFile.close();
}

void save1DVectorDoubleToBinaryFile(
    const std::vector<double>& vec,
    const std::string& filename)
{
    std::ofstream file(filename, std::ios::binary);

    if (!file)
    {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    size_t size = vec.size();
    file.write(reinterpret_cast<const char*>(&size), sizeof(size));

    file.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(double));

    file.close();
}

void save1DVectorIntToBinaryFile(
    const std::vector<int>& vec,
    const std::string& filename)
{
    std::ofstream outfile(filename, std::ios::binary);

    if (!outfile)
    {
        std::cerr << "Error opening file for writing: " << filename << std::endl;
        return;
    }

    size_t size = vec.size();
    outfile.write(reinterpret_cast<const char*>(&size), sizeof(size));

    outfile.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(int));

    outfile.close();
}

double van_der_corput(const int n, const int base)
{
    double result = 0.0;
    double f = 1.0 / base;
    int i = n;

    while (i > 0)
    {
        result += f * (i % base);
        i = i / base;
        f = f / base;
    }

    return result;
}

std::vector<std::vector<double>> generate_halton_seq(const int dim, const int n, const int skip)
{
    std::vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29}; 
    std::vector<std::vector<double>> sequence(n, std::vector<double>(dim));

    for (int i = 0; i < n; ++i)
    {
        for (int d = 0; d < dim; ++d)
        {
            sequence[i][d] = van_der_corput(i + skip, primes[d]); 
        }
    }

    return sequence;
}

std::vector<std::vector<int>> generate_halton_seq_2d_integer(const int num_points, const int x_lower_bound, const int x_upper_bound, const int y_lower_bound, const int y_upper_bound, const double rand_num)
{

    assert(rand_num >= 0.0 && rand_num <= 1.0);

    std::vector<std::vector<int>> sequence(num_points, std::vector<int>(2));

    const int skip = rand_num * 9000 + 1000;    
    for (int i = 0; i < num_points; ++i)
    {
        sequence[i][0] = std::ceil(x_lower_bound + van_der_corput(i + skip, 2) * (x_upper_bound - x_lower_bound));
        sequence[i][1] = std::ceil(y_lower_bound + van_der_corput(i + skip, 3) * (y_upper_bound - y_lower_bound));
    }

    return sequence;
}

void read_coordinate_from_file(const std::string& filename, std::vector<std::vector<double>>& x, int num_layer, int num_width)
{
    std::ifstream infile(filename);
    if (!infile) { std::cerr << "Error opening file " << filename << " in read_coordinate_from_file()" << std::endl; exit(1); }

    std::vector<std::array<double, 2>> data_coordinate;
    double value_x, value_y;
    while (infile >> value_x >> value_y)
    {
        data_coordinate.push_back({value_x, value_y});
    }
    infile.close();

    const int dof = num_layer * num_width;
    if (data_coordinate.size() != static_cast<size_t>(dof))
    {
        std::cerr << "Error: data size does not match dof" << std::endl; exit(1);
    }

    for (int j = 0; j < num_layer; ++j)
    {
        for (int i = 0; i < num_width; ++i) 
        {
            x[i + j * num_width][0] = data_coordinate[i + j * num_width][0];
            x[i + j * num_width][1] = data_coordinate[i + j * num_width][1];
        }
    }
}