#include "energy.h"
#include "config.h"
#include "constitutive.h"
#include "gradient.h"
#include "utils.h"

// example input:
// 
// ./main 20 80 10 5.0 40 1e-2 100000 1e-1 1e-1 5 0.0 2 2 2 0 0.0 1e-3 0 0.5 100 1 1234
// 
// num_layer                        20
// num_width                        80
// max_node_num                     10
// total_stretch                    5.0
// num_sampling_stretch             40
// error_tol                        1e-2
// max_time_step                    100000
// dt_no_entangle                   1e-1
// dt_entangle                      1e-1
// num_decay_dt                     5
// rand_loc                         0.0
// unit_size_x                      2
// unit_size_y                      2
// num_cross_linker_in_unit         2
// is_random_cross_linker           0
// prob_cross_linker                0.0
// error_tol_energy_stretch         1e-3
// is_random_entangle_orientation   0
// prob_entangle_orientation        0.5
// print_time_interval              100
// is_fracture                      1
// rand_seed                        1234

int main(int argc, char* argv[])
{    
    int rand_seed = 1234;            

    const std::clock_t time_start = std::clock();     
    std::clock_t time_now = std::clock();             
    double elapsed_time = static_cast<double>(time_now - time_start) / CLOCKS_PER_SEC;

    int num_layer = 16;
    int num_width = 16;

    int max_chain_num = 4;
    int max_node_num = 100;

    double total_stretch = 5.0;    

    int num_sampling_stretch = 10;
    double error_tol = 1e-5;
    int max_time_step = 1e4;
    double dt_no_entangle = 1e-1;
    double dt_entangle = 1e-2;
    int num_decay_dt = 3;
    double rand_loc = 0.05;

    int unit_size_x = 4;                
    int unit_size_y = 4;                
    int num_cross_linker_in_unit = 8;   

    int is_random_cross_linker = 0;
    double prob_cross_linker = 0.5;     

    int is_random_entangle_orientation = 0;    
    double prob_entangle_orientation = 0.5;    

    double error_tol_energy_stretch = 1e-3;    

    int is_fracture = 1;    
    int print_time_interval = 10;    

    const int min_gap_x = 2;
    const int min_gap_y = 2;

    const double stretch_crack = 5.0;                       
    const double stretch_crack_error_theshold = 0.01;        

    const double stretch_crack_init_min = 1.0;
    const double stretch_crack_init_max = 10.0;

    if (argc != 23)
    {
        if (argc < 23)
        {
            std::cout << "Too few arguments" << std::endl;
        }
        else
        {
            std::cout << "Too many arguments" << std::endl;
        }

        std::cout << "Usage: " << argv[0] << " num_layer num_width max_node_num total_stretch num_sampling_stretch error_tol max_time_step dt_no_entangle dt_entangle num_decay_dt rand_loc unit_size_x unit_size_y num_cross_linker_in_unit is_random_cross_linker prob_cross_linker is_random_entangle_orientation prob_entangle_orientation rand_seed" << std::endl; 
        return 1;
    }

    num_layer = std::atoi(argv[1]);
    num_width = std::atoi(argv[2]);
    max_node_num = std::atoi(argv[3]);
    total_stretch = std::atof(argv[4]);
    num_sampling_stretch = std::atoi(argv[5]);
    error_tol = std::atof(argv[6]);
    max_time_step = std::atoi(argv[7]);
    dt_no_entangle = std::atof(argv[8]);
    dt_entangle = std::atof(argv[9]);
    num_decay_dt = std::atoi(argv[10]);
    rand_loc = std::atof(argv[11]);
    unit_size_x = std::atoi(argv[12]);
    unit_size_y = std::atoi(argv[13]);
    num_cross_linker_in_unit = std::atoi(argv[14]);
    is_random_cross_linker = std::atoi(argv[15]);
    prob_cross_linker = std::atof(argv[16]);
    error_tol_energy_stretch = std::atof(argv[17]);
    is_random_entangle_orientation = std::atoi(argv[18]);
    prob_entangle_orientation = std::atof(argv[19]);
    print_time_interval = std::atoi(argv[20]);
    is_fracture = std::atoi(argv[21]);
    rand_seed = std::atoi(argv[22]);

    assert((max_time_step >= print_time_interval) && "max_time_step should be greater than print_time_interval");
    assert((max_time_step % print_time_interval == 0) && "max_time_step should be divisible by print_time_interval");

    std::mt19937 gen(rand_seed);     
    std::uniform_real_distribution<> distr(-1.0, 1.0);   

    int dof = num_layer * num_width;
    double total_displacement = (num_layer - 1) * (total_stretch - 1);

    std::string file_name_suffix =
        "_lay_" + std::to_string(static_cast<int>(num_layer))
        + "_wid_" + std::to_string(static_cast<int>(num_width))
        + "_maxnode_" + std::to_string(static_cast<int>(max_node_num))
        + "_stret_" + std::to_string(static_cast<int>(total_stretch))
        + "_samp_" + std::to_string(static_cast<int>(num_sampling_stretch))
        + "_tol_" + std::to_string(static_cast<int>(-log10(error_tol)))
        + "_maxstp_" + std::to_string(static_cast<int>(max_time_step))
        + "_dtne_" + std::to_string(static_cast<int>(-log10(dt_no_entangle)))
        + "_dte_" + std::to_string(static_cast<int>(-log10(dt_entangle)))
        + "_ntdcy_" + std::to_string(static_cast<int>(num_decay_dt))
        + "_randloc_" + std::to_string(static_cast<int>(rand_loc * 100))
        + "_unitx_" + std::to_string(static_cast<int>(unit_size_x))
        + "_unity_" + std::to_string(static_cast<int>(unit_size_y))
        + "_ncross_" + std::to_string(static_cast<int>(num_cross_linker_in_unit))
        + "_randcross_" + std::to_string(static_cast<int>(is_random_cross_linker))
        + "_pcross_" + std::to_string(static_cast<int>(prob_cross_linker * 100))
        + "_energytol_" + std::to_string(static_cast<int>(-log10(error_tol_energy_stretch)))
        + "_randent_" + std::to_string(static_cast<int>(is_random_entangle_orientation))
        + "_pent_" + std::to_string(static_cast<int>(prob_entangle_orientation * 100))
        + "_print_" + std::to_string(static_cast<int>(print_time_interval))
        + "_frac_" + std::to_string(static_cast<int>(is_fracture))
        + "_seed_" + std::to_string(static_cast<int>(rand_seed))
        + "_cpp";

    std::string output_file_name = "log" + file_name_suffix + ".txt";
    std::ofstream file_output(output_file_name);

    file_output << "number of layers: " << num_layer << std::endl;
    file_output << "number of width: " << num_width << std::endl;
    file_output << "total stretch: " << total_stretch << std::endl;
    file_output << "total displacement: " << total_displacement << std::endl;
    file_output << "number of sampling stretch: " << num_sampling_stretch << std::endl;
    file_output << "error tolerance: " << error_tol << std::endl;
    file_output << "max time step: " << max_time_step << std::endl;
    file_output << "time step without entangle: " << dt_no_entangle << std::endl;
    file_output << "time step with entangle: " << dt_entangle << std::endl;
    file_output << "number of decay time step: " << num_decay_dt << std::endl;
    file_output << "random location: " << rand_loc << std::endl;
    file_output << "unit size x: " << unit_size_x << std::endl;
    file_output << "unit size y: " << unit_size_y << std::endl;
    file_output << "number of cross linker in unit: " << num_cross_linker_in_unit << std::endl;
    file_output << "is random cross linker: " << is_random_cross_linker << std::endl;
    file_output << "probability of cross linker: " << prob_cross_linker << std::endl;
    file_output << "error tolerance for energy and stretch: " << error_tol_energy_stretch << std::endl;
    file_output << "is random entangle orientation: " << is_random_entangle_orientation << std::endl;
    file_output << "probability of entangle orientation: " << prob_entangle_orientation << std::endl;
    file_output << "print time interval: " << print_time_interval << std::endl;
    file_output << "is fracture: " << is_fracture << std::endl;
    file_output << "random seed: " << rand_seed << std::endl;

    std::vector<std::vector<double>> x(dof, std::vector<double>(2, 0.0));
    for (int j = 0; j < num_layer; ++j)
    {
        for (int i = 0; i < num_width; ++i)
        {
            x[i + j * num_width][0] = i;
            x[i + j * num_width][1] = j;
        }
    }

    std::vector<std::vector<int>> neighbour_matrix(dof, std::vector<int>(4, 0));
    generate_network_neighbour_matrix(num_layer, num_width, x, neighbour_matrix, is_fracture);

    std::vector<int> is_cross_linker(dof, 0);

    if (is_random_cross_linker == 0)
    {    
        initial_cross_linker_periodic_unit(num_layer, num_width, unit_size_x, unit_size_y, num_cross_linker_in_unit, is_fracture, is_cross_linker);
    }

    else if (is_random_cross_linker == 1)
    {
        initial_cross_linker_random(num_layer, num_width, prob_cross_linker, rand_seed, is_fracture, is_cross_linker);
    }

    else if (is_random_cross_linker == 2)
    {
        initial_cross_linker_quasi_random(num_layer, num_width, prob_cross_linker, rand_seed, is_cross_linker);
    }

    else if (is_random_cross_linker == 3)
    {
        initial_cross_linker_random_min_gap(num_layer, num_width, prob_cross_linker, rand_seed, min_gap_x, min_gap_y, is_fracture,is_cross_linker);
    }

    else if (is_random_cross_linker == 4)
    {
        initial_cross_linker_random_min_gap_2(num_layer, num_width, prob_cross_linker, rand_seed, is_fracture,is_cross_linker);
    }    
    else
    {
        std::cerr << "Error: unknown type of cross linker distribution" << std::endl; exit(1);
    }

    double cross_linker_ratio = compute_cross_linker_ratio_inner_domain(num_layer, num_width, is_cross_linker);
    file_output << "cross linker ratio (inner domain): " << cross_linker_ratio << std::endl;

    std::vector<int> entangle_orientation(dof, 0);

    if (is_random_entangle_orientation == 1)
    {
        for (int i = 0; i < dof; ++i)
        {
            if (distr(gen) <= 0.0)
            {
                entangle_orientation[i] = 1;
            }        
        }
    }

    std::vector<int> chain_num(dof, 0);

    std::vector<std::vector<std::vector<int>>> chain_info(dof, std::vector<std::vector<int>>(max_chain_num, std::vector<int>(max_node_num, 0)));
    int real_max_node_num = 0;
    generate_entangle_network_chain_info(num_layer, num_width, max_chain_num, max_node_num, is_cross_linker, entangle_orientation, neighbour_matrix, real_max_node_num, chain_num, chain_info);
    file_output << "max node number in one chain: " << real_max_node_num << std::endl;
    if (real_max_node_num > max_node_num)
    {
        std::cerr << "Error: real_max_node_num is greater than max_node_num" << std::endl; exit(1);
    }

    for (int j = 0; j < num_layer; ++j)
    {
        for (int i = 0; i < num_width; ++i)
        {

            bool is_b = (j == 0);
            bool is_t = (j == (num_layer - 1));
            bool is_l = (i == 0);
            bool is_r = (i == (num_width - 1));

            if (!is_b && !is_t && !is_l && !is_r)
            {
                x[i + j * num_width][0] += distr(gen) * rand_loc;
                x[i + j * num_width][1] += distr(gen) * rand_loc;                
            }
        }
    }

    const std::vector<std::vector<double>> x_init = x;

    std::vector<std::vector<double>> chain_init_len(dof, std::vector<double>(max_chain_num, 0.0));
    std::vector<std::vector<double>> dist_matrix_init;
    compute_chain_init_len(num_layer, num_width, x, neighbour_matrix, chain_num, chain_info, dist_matrix_init, chain_init_len);

if (is_fracture == 0)
{
    std::vector<double> stretch_list(num_sampling_stretch + 1);
    for (int i = 0; i < num_sampling_stretch + 1; ++i)
    {
        stretch_list[i] = 1 + i * (total_stretch - 1) / num_sampling_stretch;
    }
    stretch_list.erase(stretch_list.begin());

    std::vector<double> displacement_list(num_sampling_stretch);
    for (int i = 0; i < num_sampling_stretch; ++i)
    {
        displacement_list[i] = (num_layer - 1) * (stretch_list[i] - 1);

    }

    std::vector<std::vector<std::vector<double>>> x_entangle_list(num_sampling_stretch, std::vector<std::vector<double>>(dof, std::vector<double>(2, 0.0)));
    std::vector<std::vector<std::vector<double>>> x_no_entangle_list(num_sampling_stretch, std::vector<std::vector<double>>(dof, std::vector<double>(2, 0.0)));

    std::vector<double> energy_entangle_list(num_sampling_stretch, 0.0);
    std::vector<double> energy_no_entangle_list(num_sampling_stretch, 0.0);

    std::vector<double> max_stretch_entangle_list(num_sampling_stretch, 0.0);
    std::vector<double> max_stretch_no_entangle_list(num_sampling_stretch, 0.0);

    for (int i_num_sampling = 0; i_num_sampling < num_sampling_stretch; ++i_num_sampling)
    {
        file_output << "\n";
        file_output << "---------- \n";
        file_output << "start solve network at number of sampling " << i_num_sampling + 1 << " / " << num_sampling_stretch << std::endl;
        file_output << "---------- \n";
        file_output << "STEP 1: solve network without entangle \n";
        file_output << "the result will be taken as initial guess to solve network with entangle \n";
        file_output << "---------- \n";

        file_output << "displacement " << displacement_list[i_num_sampling] << std::endl;
        file_output << "stretch " << stretch_list[i_num_sampling] << std::endl;

        std::vector<std::vector<double>> xb = x_init;
        for (int i = 0; i < num_width; ++i)
        {
            xb[dof - i - 1][1] += displacement_list[i_num_sampling];
        }

        double y_upper_boundary = x_init[dof-1][1] + displacement_list[i_num_sampling];
        for (int j = 0; j < num_layer; ++j)
        {
            for (int i = 0; i < num_width; ++i)
            {
                x[i + j * num_width][0] = x_init[i + j * num_width][0]; 
                x[i + j * num_width][1] = y_upper_boundary * j / (num_layer - 1.0);
            }
        }

        time_now = std::clock();
        x_no_entangle_list[i_num_sampling] = solve_optimization_elasticity(x, xb, false, dt_no_entangle, num_decay_dt, max_time_step, print_time_interval, num_width, num_layer, error_tol, neighbour_matrix, dist_matrix_init, chain_num, chain_init_len, chain_info, file_output);
        x = x_no_entangle_list[i_num_sampling];

        elapsed_time = static_cast<double>(std::clock() - time_now) / CLOCKS_PER_SEC;
        file_output << std::scientific << std::setprecision(2) << "elapsed time in solving network without entangle " << elapsed_time << " sec" << std::endl;

        double total_energy = compute_total_energy_no_entangle_inner_domain(x, neighbour_matrix, num_layer, num_width, dist_matrix_init);
        energy_no_entangle_list[i_num_sampling] = total_energy;
        double max_stretch = compute_max_stretch_no_entangle_inner_domain(x, neighbour_matrix, num_layer, num_width, dist_matrix_init);
        max_stretch_no_entangle_list[i_num_sampling] = max_stretch;
        file_output << "total energy without entangle " << total_energy << std::endl;
        file_output << "max stretch without entangle " << max_stretch << std::endl;
        file_output << "---------- \n";

        file_output << "STEP 2: solve network with entangle \n";
        file_output << "---------- \n";

        time_now = std::clock();
        x_entangle_list[i_num_sampling] = solve_optimization_elasticity(x, xb, true, dt_entangle, num_decay_dt, max_time_step, print_time_interval, num_width, num_layer, error_tol, neighbour_matrix, dist_matrix_init, chain_num, chain_init_len, chain_info, file_output);
        x = x_entangle_list[i_num_sampling];

        elapsed_time = static_cast<double>(std::clock() - time_now) / CLOCKS_PER_SEC;
        file_output << std::scientific << std::setprecision(2) << "elapsed time in solving network with entangle " << elapsed_time << " sec" << std::endl;

        total_energy = compute_total_energy_entangle_inner_domain(x, neighbour_matrix, chain_num, chain_init_len, chain_info, num_layer, num_width);
        energy_entangle_list[i_num_sampling] = total_energy;
        max_stretch = compute_max_stretch_entangle_inner_domain(x, neighbour_matrix, chain_num, chain_init_len, chain_info, num_layer, num_width);
        max_stretch_entangle_list[i_num_sampling] = max_stretch;
        file_output << "total energy with entangle " << total_energy << std::endl;
        file_output << "max stretch with entangle " << max_stretch << std::endl;
        file_output << "---------- \n";
    }

    elapsed_time = static_cast<double>(std::clock() - time_start) / CLOCKS_PER_SEC;
    file_output << "---------- \n";
    file_output << "computation completed" << std::endl;
    file_output << "total elapsed time " << elapsed_time << " sec" << std::endl;
    file_output << "---------- \n";    

    file_output << "---------- \n";
    file_output << "start writing variables to binary file" << std::endl;

    save2DVectorDoubleToBinaryFile(x_init, "x_init" + file_name_suffix + ".bin");

    save3DVectorDoubleToBinaryFile(x_entangle_list, "x_entangle" + file_name_suffix + ".bin");

    save3DVectorDoubleToBinaryFile(x_no_entangle_list, "x_no_entangle" + file_name_suffix + ".bin");

    save1DVectorDoubleToBinaryFile(energy_entangle_list, "energy_entangle" + file_name_suffix + ".bin");

    save1DVectorDoubleToBinaryFile(energy_no_entangle_list, "energy_no_entangle" + file_name_suffix + ".bin");

    save1DVectorDoubleToBinaryFile(max_stretch_entangle_list, "max_stretch_entangle" + file_name_suffix + ".bin");

    save1DVectorDoubleToBinaryFile(max_stretch_no_entangle_list, "max_stretch_no_entangle" + file_name_suffix + ".bin");

    save1DVectorIntToBinaryFile(is_cross_linker, "is_cross_linker" + file_name_suffix + ".bin");

    save1DVectorIntToBinaryFile(entangle_orientation, "entangle_orientation" + file_name_suffix + ".bin");

    save1DVectorIntToBinaryFile(chain_num, "chain_num" + file_name_suffix + ".bin");
    save2DVectorDoubleToBinaryFile(chain_init_len, "chain_init_len" + file_name_suffix + ".bin");

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < max_chain_num; ++j)
        {
            for (int k = 0; k < max_node_num; ++k)
            {
                if (k >= int(chain_info[i][j].size()))
                {
                    chain_info[i][j].push_back(-1);
                }
                chain_info[i][j][k] += 1;
            }
        }
    }
    save3DVectorIntToBinaryFile(chain_info, "chain_info" + file_name_suffix + ".bin");

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            neighbour_matrix[i][j] += 1;
        }
    }
    save2DVectorIntToBinaryFile(neighbour_matrix, "neighbour_matrix" + file_name_suffix + ".bin");

    file_output << "variables are saved " << std::endl;
    file_output << "---------- \n";
    file_output.close();

    if (num_sampling_stretch == 1)
    {
        std::string result_file_name_single_sampling = "result_single_sampling" + file_name_suffix + ".txt";
        std::ofstream result_output(result_file_name_single_sampling);

        result_output << std::scientific << std::setprecision(6) << max_stretch_no_entangle_list[0] << std::endl;
        result_output << std::scientific << std::setprecision(6) << energy_no_entangle_list[0] << std::endl;

        result_output << std::scientific << std::setprecision(6) << max_stretch_entangle_list[0] << std::endl;
        result_output << std::scientific << std::setprecision(6) << energy_entangle_list[0] << std::endl;

        result_output.close();
    }

    return 0;
}

else if (is_fracture == 1)
{
    std::vector<std::vector<double>> x_no_entangle_crack = x_init;
    std::vector<std::vector<double>> x_entangle_crack = x_init;

    std::vector<double> stretch_list_no_entangle;
    std::vector<double> max_stretch_list_no_entangle;

    file_output << "\n";
    file_output << "-------------------------------------------------- \n";
    file_output << "start solving fracture network with no entangle \n";
    file_output << "-------------------------------------------------- \n";
    file_output << "\n";
    {        

        double stretch_c_min = stretch_crack_init_min;
        double stretch_c_max = stretch_crack_init_max;
        bool is_bisection_converge = false;
        int num_iter_bisection = 0;        
        while (!is_bisection_converge)
        {
            double stretch_c = 0.5 * (stretch_c_min + stretch_c_max);
            stretch_list_no_entangle.push_back(stretch_c);

            file_output << "\n";
            file_output << "---------- \n";
            file_output << "find crack stretch for network without entangle in number of iteration: " << num_iter_bisection + 1 << " \n";
            file_output << std::scientific << std::setprecision(6) << "try stretch: " << stretch_c << std::endl;
            file_output << "---------- \n";

            double disp_c = (num_layer - 1) * (stretch_c - 1);

            std::vector<std::vector<double>> xb = x_init;
            for (int i = 0; i < num_width; ++i)
            {
                xb[dof - i - 1][1] += disp_c;
            }

            double y_upper_boundary = x_init[dof-1][1] + disp_c;
            for (int j = 0; j < num_layer; ++j)
            {
                for (int i = 0; i < num_width; ++i)
                {
                    x[i + j * num_width][0] = x_init[i + j * num_width][0]; 
                    x[i + j * num_width][1] = y_upper_boundary * j / (num_layer - 1.0);
                }
            }

            time_now = std::clock();

            x = solve_optimization_fracture(x, xb, false, dt_no_entangle, num_decay_dt, max_time_step, print_time_interval, num_width, num_layer, error_tol, neighbour_matrix, dist_matrix_init, chain_num, chain_init_len, chain_info, file_output, 1, stretch_crack);

            elapsed_time = static_cast<double>(std::clock() - time_now) / CLOCKS_PER_SEC;
            file_output << "---------- \n";
            file_output << "elapsed time in solving network without entangle " << elapsed_time << " sec" << std::endl;

            double max_stretch = compute_max_stretch_no_entangle_inner_domain(x, neighbour_matrix, num_layer, num_width, dist_matrix_init);
            max_stretch_list_no_entangle.push_back(max_stretch);

            if ((std::abs(max_stretch - stretch_crack) <= stretch_crack_error_theshold) || (std::abs(stretch_c_max - stretch_c_min) <= stretch_crack_error_theshold))
            {
                is_bisection_converge = true;
                x_no_entangle_crack = x;

                file_output << std::scientific << std::setprecision(6) << "max stretch " << max_stretch << " is close to the stretch for crack " << stretch_crack << "; network without entangle is solved. \n";
                file_output << "stretch for crack: " << stretch_c << std::endl;
                file_output << "---------- \n";

                double elapsed_time = static_cast<double>(std::clock() - time_start) / CLOCKS_PER_SEC;
                file_output << "---------- \n";
                file_output << "end solving fracture network without entangle \n";
                file_output << "total elapsed time " << elapsed_time << " sec" << std::endl;
                file_output << "---------- \n";    
            }
            else if (max_stretch < stretch_crack)
            {
                file_output << "---------- \n";
                file_output << std::scientific << std::setprecision(6) << "max stretch " << max_stretch << " is less than the stretch for crack " << stretch_crack << "; increase the stretch for crack. \n";
                file_output << "---------- \n";
                stretch_c_min = stretch_c;
            }
            else if (max_stretch > stretch_crack)
            {
                file_output << "---------- \n";
                file_output << std::scientific << std::setprecision(6) << "max stretch " << max_stretch << " is greater than the stretch for crack " << stretch_crack << "; decrease the stretch for crack. \n";
                file_output << "---------- \n";
                stretch_c_max = stretch_c;
            }

            num_iter_bisection += 1;
            if (num_iter_bisection > 20)
            {
                file_output << "---------- \n";
                file_output << "bisection method does not converge; stop the computation. \n";
                file_output << "---------- \n";
                return 0;
            }
        }
    }    

    std::vector<double> stretch_list_entangle;
    std::vector<double> max_stretch_list_entangle;

    file_output << "\n";
    file_output << "-------------------------------------------------- \n";
    file_output << "start solving fracture network with entangle \n";
    file_output << "-------------------------------------------------- \n";
    file_output << "\n";    
    {

        double stretch_c_min = stretch_crack_init_min;
        double stretch_c_max = stretch_crack_init_max;
        bool is_bisection_converge = false;
        int num_iter_bisection = 0;
        while (!is_bisection_converge)
        {
            double stretch_c = 0.5 * (stretch_c_min + stretch_c_max);
            stretch_list_entangle.push_back(stretch_c);

            file_output << "\n";
            file_output << "---------- \n";
            file_output << "find crack stretch for network with entangle in number of iteration: " << num_iter_bisection + 1 << " \n";
            file_output << std::scientific << std::setprecision(6) << "try stretch: " << stretch_c << std::endl;
            file_output << "---------- \n";

            double disp_c = (num_layer - 1) * (stretch_c - 1);

            std::vector<std::vector<double>> xb = x_init;
            for (int i = 0; i < num_width; ++i)
            {
                xb[dof - i - 1][1] += disp_c;
            }

            double y_upper_boundary = x_init[dof-1][1] + disp_c;
            for (int j = 0; j < num_layer; ++j)
            {
                for (int i = 0; i < num_width; ++i)
                {
                    x[i + j * num_width][0] = x_init[i + j * num_width][0]; 
                    x[i + j * num_width][1] = y_upper_boundary * j / (num_layer - 1.0);
                }
            }

            time_now = std::clock();

            x = solve_optimization_fracture(x, xb, true, dt_no_entangle, num_decay_dt, max_time_step, print_time_interval, num_width, num_layer, error_tol, neighbour_matrix, dist_matrix_init, chain_num, chain_init_len, chain_info, file_output, 1, stretch_crack);

            elapsed_time = static_cast<double>(std::clock() - time_now) / CLOCKS_PER_SEC;
            file_output << "---------- \n";
            file_output << "elapsed time in solving network with entangle " << elapsed_time << " sec" << std::endl;

            double max_stretch = compute_max_stretch_entangle_inner_domain(x, neighbour_matrix, chain_num, chain_init_len, chain_info, num_layer, num_width);
            max_stretch_list_entangle.push_back(max_stretch);

            if ((std::abs(max_stretch - stretch_crack) <= stretch_crack_error_theshold) || (std::abs(stretch_c_max - stretch_c_min) <= stretch_crack_error_theshold))
            {
                is_bisection_converge = true;
                x_entangle_crack = x;

                file_output << std::scientific << std::setprecision(6) << "max stretch " << max_stretch << " is close to the stretch for crack " << stretch_crack << "; network with entangle is solved. \n";
                file_output << "stretch for crack: " << stretch_c << std::endl;
                file_output << "---------- \n";

                double elapsed_time = static_cast<double>(std::clock() - time_start) / CLOCKS_PER_SEC;
                file_output << "---------- \n";
                file_output << "end solving fracture network with entangle \n";
                file_output << "total elapsed time " << elapsed_time << " sec" << std::endl;
                file_output << "---------- \n";
            }
            else if (max_stretch < stretch_crack)
            {
                file_output << "---------- \n";
                file_output << std::scientific << std::setprecision(6) << "max stretch " << max_stretch << " is less than the stretch for crack " << stretch_crack << "; increase the stretch for crack. \n";
                file_output << "---------- \n";
                stretch_c_min = stretch_c;
            }
            else if (max_stretch > stretch_crack)
            {
                file_output << "---------- \n";
                file_output << std::scientific << std::setprecision(6) << "max stretch " << max_stretch << " is greater than the stretch for crack " << stretch_crack << "; decrease the stretch for crack. \n";
                file_output << "---------- \n";
                stretch_c_max = stretch_c;
            }

            num_iter_bisection += 1;
            if (num_iter_bisection > 20)
            {
                file_output << "---------- \n";
                file_output << "bisection method does not converge; stop the computation. \n";
                file_output << "---------- \n";
                return 0;
            }
        }
    }

    file_output << "---------- \n";
    file_output << "start writing variables to binary file" << std::endl;

    save2DVectorDoubleToBinaryFile(x_init, "x_init" + file_name_suffix + ".bin");

    save2DVectorDoubleToBinaryFile(x_entangle_crack, "x_entangle_crack" + file_name_suffix + ".bin");

    save2DVectorDoubleToBinaryFile(x_no_entangle_crack, "x_no_entangle_crack" + file_name_suffix + ".bin");

    save1DVectorIntToBinaryFile(is_cross_linker, "is_cross_linker" + file_name_suffix + ".bin");

    save1DVectorIntToBinaryFile(chain_num, "chain_num" + file_name_suffix + ".bin");
    save2DVectorDoubleToBinaryFile(chain_init_len, "chain_init_len" + file_name_suffix + ".bin");

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < max_chain_num; ++j)
        {
            for (int k = 0; k < max_node_num; ++k)
            {
                if (k >= int(chain_info[i][j].size()))
                {
                    chain_info[i][j].push_back(-1);
                }
                chain_info[i][j][k] += 1;
            }
        }
    }
    save3DVectorIntToBinaryFile(chain_info, "chain_info" + file_name_suffix + ".bin");

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            neighbour_matrix[i][j] += 1;
        }
    }
    save2DVectorIntToBinaryFile(neighbour_matrix, "neighbour_matrix" + file_name_suffix + ".bin");

    save1DVectorDoubleToBinaryFile(stretch_list_no_entangle, "stretch_list_no_entangle" + file_name_suffix + ".bin");
    save1DVectorDoubleToBinaryFile(max_stretch_list_no_entangle, "max_stretch_list_no_entangle" + file_name_suffix + ".bin");

    save1DVectorDoubleToBinaryFile(stretch_list_entangle, "stretch_list_entangle" + file_name_suffix + ".bin");
    save1DVectorDoubleToBinaryFile(max_stretch_list_entangle, "max_stretch_list_entangle" + file_name_suffix + ".bin");

    file_output << "variables are saved " << std::endl;
    file_output << "---------- \n";
    file_output.close();

    return 0;
}

else if (is_fracture == 2)
{
    std::vector<double> stretch_list(num_sampling_stretch + 1);
    for (int i = 0; i < num_sampling_stretch + 1; ++i)
    {
        stretch_list[i] = 1 + i * (total_stretch - 1) / num_sampling_stretch;
    }
    stretch_list.erase(stretch_list.begin());

    std::vector<double> displacement_list(num_sampling_stretch);
    for (int i = 0; i < num_sampling_stretch; ++i)
    {
        displacement_list[i] = (num_layer - 1) * (stretch_list[i] - 1);
    }

    std::vector<std::vector<std::vector<double>>> x_entangle_list(num_sampling_stretch, std::vector<std::vector<double>>(dof, std::vector<double>(2, 0.0)));
    std::vector<std::vector<std::vector<double>>> x_no_entangle_list(num_sampling_stretch, std::vector<std::vector<double>>(dof, std::vector<double>(2, 0.0)));

    std::vector<double> energy_entangle_list(num_sampling_stretch, 0.0);
    std::vector<double> energy_no_entangle_list(num_sampling_stretch, 0.0);

    std::vector<double> max_stretch_entangle_list(num_sampling_stretch, 0.0);
    std::vector<double> max_stretch_no_entangle_list(num_sampling_stretch, 0.0);

    for (int i_num_sampling = 0; i_num_sampling < num_sampling_stretch; ++i_num_sampling)
    {
        file_output << "\n";
        file_output << "---------- \n";
        file_output << "start solve network at number of sampling " << i_num_sampling + 1 << " / " << num_sampling_stretch << std::endl;
        file_output << "---------- \n";
        file_output << "STEP 1: solve network without entangle \n";
        file_output << "the result will be taken as initial guess to solve network with entangle \n";
        file_output << "---------- \n";

        file_output << "displacement " << displacement_list[i_num_sampling] << std::endl;
        file_output << "stretch " << stretch_list[i_num_sampling] << std::endl;

        std::vector<std::vector<double>> xb = x_init;
        for (int i = 0; i < num_width; ++i)
        {
            xb[dof - i - 1][1] += displacement_list[i_num_sampling];
        }

        double y_upper_boundary = x_init[dof-1][1] + displacement_list[i_num_sampling];
        for (int j = 0; j < num_layer; ++j)
        {
            for (int i = 0; i < num_width; ++i)
            {
                x[i + j * num_width][0] = x_init[i + j * num_width][0]; 
                x[i + j * num_width][1] = y_upper_boundary * j / (num_layer - 1.0);
            }
        }

        time_now = std::clock();
        x = solve_optimization_fracture(x, xb, false, dt_no_entangle, num_decay_dt, max_time_step, print_time_interval, num_width, num_layer, error_tol, neighbour_matrix, dist_matrix_init, chain_num, chain_init_len, chain_info, file_output, 0, stretch_crack);
        x_no_entangle_list[i_num_sampling] = x;

        elapsed_time = static_cast<double>(std::clock() - time_now) / CLOCKS_PER_SEC;
        file_output << std::scientific << std::setprecision(2) << "elapsed time in solving network without entangle " << elapsed_time << " sec" << std::endl;

        double total_energy = compute_total_energy_no_entangle_inner_domain(x, neighbour_matrix, num_layer, num_width, dist_matrix_init);
        energy_no_entangle_list[i_num_sampling] = total_energy;
        double max_stretch = compute_max_stretch_no_entangle_inner_domain(x, neighbour_matrix, num_layer, num_width, dist_matrix_init);
        max_stretch_no_entangle_list[i_num_sampling] = max_stretch;
        file_output << "total energy without entangle " << total_energy << std::endl;
        file_output << "max stretch without entangle " << max_stretch << std::endl;
        file_output << "---------- \n";

        file_output << "STEP 2: solve network with entangle \n";
        file_output << "---------- \n";

        xb = x_init;
        for (int i = 0; i < num_width; ++i)
        {
            xb[dof - i - 1][1] += displacement_list[i_num_sampling];
        }

        y_upper_boundary = x_init[dof-1][1] + displacement_list[i_num_sampling];
        for (int j = 0; j < num_layer; ++j)
        {
            for (int i = 0; i < num_width; ++i)
            {
                x[i + j * num_width][0] = x_init[i + j * num_width][0]; 
                x[i + j * num_width][1] = y_upper_boundary * j / (num_layer - 1.0);
            }
        }

        time_now = std::clock();
        x = solve_optimization_fracture(x, xb, true, dt_no_entangle, num_decay_dt, max_time_step, print_time_interval, num_width, num_layer, error_tol, neighbour_matrix, dist_matrix_init, chain_num, chain_init_len, chain_info, file_output, 0, stretch_crack);
        x_entangle_list[i_num_sampling] = x;

        elapsed_time = static_cast<double>(std::clock() - time_now) / CLOCKS_PER_SEC;
        file_output << std::scientific << std::setprecision(2) << "elapsed time in solving network with entangle " << elapsed_time << " sec" << std::endl;

        total_energy = compute_total_energy_entangle_inner_domain(x, neighbour_matrix, chain_num, chain_init_len, chain_info, num_layer, num_width);
        energy_entangle_list[i_num_sampling] = total_energy;
        max_stretch = compute_max_stretch_entangle_inner_domain(x, neighbour_matrix, chain_num, chain_init_len, chain_info, num_layer, num_width);
        max_stretch_entangle_list[i_num_sampling] = max_stretch;
        file_output << "total energy with entangle " << total_energy << std::endl;
        file_output << "max stretch with entangle " << max_stretch << std::endl;
        file_output << "---------- \n";
    }

    elapsed_time = static_cast<double>(std::clock() - time_start) / CLOCKS_PER_SEC;
    file_output << "---------- \n";
    file_output << "computation completed" << std::endl;
    file_output << "total elapsed time " << elapsed_time << " sec" << std::endl;
    file_output << "---------- \n";    

    file_output << "---------- \n";
    file_output << "start writing variables to binary file" << std::endl;

    save2DVectorDoubleToBinaryFile(x_init, "x_init" + file_name_suffix + ".bin");

    save3DVectorDoubleToBinaryFile(x_entangle_list, "x_entangle" + file_name_suffix + ".bin");

    save3DVectorDoubleToBinaryFile(x_no_entangle_list, "x_no_entangle" + file_name_suffix + ".bin");

    save1DVectorDoubleToBinaryFile(energy_entangle_list, "energy_entangle" + file_name_suffix + ".bin");

    save1DVectorDoubleToBinaryFile(energy_no_entangle_list, "energy_no_entangle" + file_name_suffix + ".bin");

    save1DVectorDoubleToBinaryFile(max_stretch_entangle_list, "max_stretch_entangle" + file_name_suffix + ".bin");

    save1DVectorDoubleToBinaryFile(max_stretch_no_entangle_list, "max_stretch_no_entangle" + file_name_suffix + ".bin");

    save1DVectorIntToBinaryFile(is_cross_linker, "is_cross_linker" + file_name_suffix + ".bin");

    save1DVectorIntToBinaryFile(entangle_orientation, "entangle_orientation" + file_name_suffix + ".bin");

    save1DVectorIntToBinaryFile(chain_num, "chain_num" + file_name_suffix + ".bin");
    save2DVectorDoubleToBinaryFile(chain_init_len, "chain_init_len" + file_name_suffix + ".bin");

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < max_chain_num; ++j)
        {
            for (int k = 0; k < max_node_num; ++k)
            {
                if (k >= int(chain_info[i][j].size()))
                {
                    chain_info[i][j].push_back(-1);
                }
                chain_info[i][j][k] += 1;
            }
        }
    }
    save3DVectorIntToBinaryFile(chain_info, "chain_info" + file_name_suffix + ".bin");

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            neighbour_matrix[i][j] += 1;
        }
    }
    save2DVectorIntToBinaryFile(neighbour_matrix, "neighbour_matrix" + file_name_suffix + ".bin");

    file_output << "variables are saved " << std::endl;
    file_output << "---------- \n";
    file_output.close();

    if (num_sampling_stretch == 1)
    {
        std::string result_file_name_single_sampling = "result_single_sampling" + file_name_suffix + ".txt";
        std::ofstream result_output(result_file_name_single_sampling);

        result_output << std::scientific << std::setprecision(6) << max_stretch_no_entangle_list[0] << std::endl;
        result_output << std::scientific << std::setprecision(6) << energy_no_entangle_list[0] << std::endl;

        result_output << std::scientific << std::setprecision(6) << max_stretch_entangle_list[0] << std::endl;
        result_output << std::scientific << std::setprecision(6) << energy_entangle_list[0] << std::endl;

        result_output.close();
    }

    return 0;
}

else
{
    std::cerr << "invalid input for is_fracture" << std::endl;
    return 1;
}

}