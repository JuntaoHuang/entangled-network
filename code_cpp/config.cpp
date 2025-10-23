#include "config.h"
#include <cassert>

void generate_entangle_network(
    int num_layer,
    int num_width,
    int max_chain_num,
    int max_node_num,
    std::vector<int>& chain_num,
    std::vector<std::vector<std::vector<int>>>& chain_info)
{

    assert(num_layer % 2 == 0);
    assert(num_width % 2 == 0);

    int dof = num_layer * num_width;
    chain_num.assign(dof, 0);
    chain_info.assign(dof, std::vector<std::vector<int>>(max_chain_num, std::vector<int>(max_node_num, -1)));

    for (int j = 0; j < num_layer; ++j)
    {
        for (int i = 0; i < num_width; ++i)
        {
            int node_index = local_to_global_index(i, j, num_width);

            int node_l_index = local_to_global_index(i - 1, j, num_width);
            int node_r_index = local_to_global_index(i + 1, j, num_width);
            int node_b_index = local_to_global_index(i, j - 1, num_width);
            int node_t_index = local_to_global_index(i, j + 1, num_width);

            int node_br_index = local_to_global_index(i + 1, j - 1, num_width);
            int node_tl_index = local_to_global_index(i - 1, j + 1, num_width);

            bool is_b = (j == 0);
            bool is_t = (j == (num_layer - 1));
            bool is_l = (i == 0);
            bool is_r = (i == (num_width - 1));

            if (is_b && is_l)
            {
                chain_num[node_index] = 2;

                chain_info[node_index][0][0] = node_index;
                chain_info[node_index][0][1] = node_t_index;

                chain_info[node_index][1][0] = node_index;
                chain_info[node_index][1][1] = node_r_index;
            }

            else if (is_b && is_r)
            {
                chain_num[node_index] = 2;

                chain_info[node_index][0][0] = node_index;
                chain_info[node_index][0][1] = node_t_index;

                chain_info[node_index][1][0] = node_index;
                chain_info[node_index][1][1] = node_l_index;
            }

            else if (is_t && is_l)
            {
                chain_num[node_index] = 2;

                chain_info[node_index][0][0] = node_index;
                chain_info[node_index][0][1] = node_b_index;

                chain_info[node_index][1][0] = node_index;
                chain_info[node_index][1][1] = node_r_index;
            }

            else if (is_t && is_r)
            {
                chain_num[node_index] = 2;

                chain_info[node_index][0][0] = node_index;
                chain_info[node_index][0][1] = node_b_index;

                chain_info[node_index][1][0] = node_index;
                chain_info[node_index][1][1] = node_l_index;
            }

            else if (is_b && !is_l && !is_r)
            {
                chain_num[node_index] = 3;

                chain_info[node_index][0][0] = node_index;
                chain_info[node_index][0][1] = node_l_index;

                chain_info[node_index][1][0] = node_index;
                chain_info[node_index][1][1] = node_r_index;

                chain_info[node_index][2][0] = node_index;
                chain_info[node_index][2][1] = node_t_index;
                if (i % 2 == 1)
                {
                    chain_info[node_index][2][2] = node_tl_index;
                }
            }

            else if (is_t && !is_l && !is_r)
            {
                chain_num[node_index] = 3;

                chain_info[node_index][0][0] = node_index;
                chain_info[node_index][0][1] = node_l_index;

                chain_info[node_index][1][0] = node_index;
                chain_info[node_index][1][1] = node_r_index;

                chain_info[node_index][2][0] = node_index;
                chain_info[node_index][2][1] = node_b_index;
                if (i % 2 == 0)
                {
                    chain_info[node_index][2][2] = node_br_index;
                }
            }

            else if (is_l && !is_b && !is_t)
            {
                chain_num[node_index] = 3;

                chain_info[node_index][0][0] = node_index;
                chain_info[node_index][0][1] = node_t_index;

                chain_info[node_index][1][0] = node_index;
                chain_info[node_index][1][1] = node_b_index;

                chain_info[node_index][2][0] = node_index;
                chain_info[node_index][2][1] = node_r_index;
                if (j % 2 == 1)
                {
                    chain_info[node_index][2][2] = node_br_index;
                }
            }

            else if (is_r && !is_b && !is_t)
            {
                chain_num[node_index] = 3;

                chain_info[node_index][0][0] = node_index;
                chain_info[node_index][0][1] = node_t_index;

                chain_info[node_index][1][0] = node_index;
                chain_info[node_index][1][1] = node_b_index;

                chain_info[node_index][2][0] = node_index;
                chain_info[node_index][2][1] = node_l_index;
                if (j % 2 == 0)
                {
                    chain_info[node_index][2][2] = node_tl_index;
                }
            }

            else if (!is_t && !is_b && !is_l && !is_r)
            {

                if ((i % 2 == 0 && j % 2 == 0) || (i % 2 == 1 && j % 2 == 1))
                {
                    chain_num[node_index] = 2;

                    chain_info[node_index][0][0] = node_l_index;
                    chain_info[node_index][0][1] = node_index;
                    chain_info[node_index][0][2] = node_b_index;

                    chain_info[node_index][1][0] = node_t_index;
                    chain_info[node_index][1][1] = node_index;
                    chain_info[node_index][1][2] = node_r_index;
                }

                else
                {
                    chain_num[node_index] = 4;

                    if (i == 1 && j == num_layer - 2)
                    {
                        chain_info[node_index][0][0] = node_index;
                        chain_info[node_index][0][1] = node_t_index;

                        chain_info[node_index][1][0] = node_index;
                        chain_info[node_index][1][1] = node_b_index;
                        chain_info[node_index][1][2] = node_br_index;

                        chain_info[node_index][2][0] = node_index;
                        chain_info[node_index][2][1] = node_l_index;

                        chain_info[node_index][3][0] = node_index;
                        chain_info[node_index][3][1] = node_r_index;
                        chain_info[node_index][3][2] = node_br_index;
                    }

                    else if (i == num_width - 2 && j == 1)
                    {
                        chain_info[node_index][0][0] = node_index;
                        chain_info[node_index][0][1] = node_t_index;
                        chain_info[node_index][0][2] = node_tl_index;

                        chain_info[node_index][1][0] = node_index;
                        chain_info[node_index][1][1] = node_l_index;
                        chain_info[node_index][1][2] = node_tl_index;

                        chain_info[node_index][2][0] = node_index;
                        chain_info[node_index][2][1] = node_b_index;

                        chain_info[node_index][3][0] = node_index;
                        chain_info[node_index][3][1] = node_r_index;
                    }

                    else if (i == 1 && j < num_layer - 2)
                    {
                        chain_info[node_index][0][0] = node_index;
                        chain_info[node_index][0][1] = node_t_index;
                        chain_info[node_index][0][2] = node_tl_index;

                        chain_info[node_index][1][0] = node_index;
                        chain_info[node_index][1][1] = node_l_index;

                        chain_info[node_index][2][0] = node_index;
                        chain_info[node_index][2][1] = node_b_index;
                        chain_info[node_index][2][2] = node_br_index;

                        chain_info[node_index][3][0] = node_index;
                        chain_info[node_index][3][1] = node_r_index;
                        chain_info[node_index][3][2] = node_br_index;
                    }

                    else if (i == num_width - 2 && j > 1)
                    {
                        chain_info[node_index][0][0] = node_index;
                        chain_info[node_index][0][1] = node_t_index;
                        chain_info[node_index][0][2] = node_tl_index;

                        chain_info[node_index][1][0] = node_index;
                        chain_info[node_index][1][1] = node_l_index;
                        chain_info[node_index][1][2] = node_tl_index;

                        chain_info[node_index][2][0] = node_index;
                        chain_info[node_index][2][1] = node_b_index;
                        chain_info[node_index][2][2] = node_br_index;

                        chain_info[node_index][3][0] = node_index;
                        chain_info[node_index][3][1] = node_r_index;
                    }

                    else if (i < num_width - 2 && j == 1)
                    {
                        chain_info[node_index][0][0] = node_index;
                        chain_info[node_index][0][1] = node_t_index;
                        chain_info[node_index][0][2] = node_tl_index;

                        chain_info[node_index][1][0] = node_index;
                        chain_info[node_index][1][1] = node_l_index;
                        chain_info[node_index][1][2] = node_tl_index;

                        chain_info[node_index][2][0] = node_index;
                        chain_info[node_index][2][1] = node_b_index;

                        chain_info[node_index][3][0] = node_index;
                        chain_info[node_index][3][1] = node_r_index;
                        chain_info[node_index][3][2] = node_br_index;
                    }

                    else if (i > 1 && j == num_layer - 2)
                    {
                        chain_info[node_index][0][0] = node_index;
                        chain_info[node_index][0][1] = node_t_index;

                        chain_info[node_index][1][0] = node_index;
                        chain_info[node_index][1][1] = node_l_index;
                        chain_info[node_index][1][2] = node_tl_index;

                        chain_info[node_index][2][0] = node_index;
                        chain_info[node_index][2][1] = node_b_index;
                        chain_info[node_index][2][2] = node_br_index;

                        chain_info[node_index][3][0] = node_index;
                        chain_info[node_index][3][1] = node_r_index;
                        chain_info[node_index][3][2] = node_br_index;
                    }

                    else
                    {
                        chain_info[node_index][0][0] = node_index;
                        chain_info[node_index][0][1] = node_t_index;
                        chain_info[node_index][0][2] = node_tl_index;

                        chain_info[node_index][1][0] = node_index;
                        chain_info[node_index][1][1] = node_l_index;
                        chain_info[node_index][1][2] = node_tl_index;

                        chain_info[node_index][2][0] = node_index;
                        chain_info[node_index][2][1] = node_b_index;
                        chain_info[node_index][2][2] = node_br_index;

                        chain_info[node_index][3][0] = node_index;
                        chain_info[node_index][3][1] = node_r_index;
                        chain_info[node_index][3][2] = node_br_index;
                    }
                }
            }
        }
    }

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < chain_num[i]; ++j)
        {
            std::vector<int> node_in_chain = chain_info[i][j];
            node_in_chain = remove_const(node_in_chain, -1);

            chain_info[i][j] = node_in_chain;
        }
    }
}

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
    std::vector<std::vector<std::vector<int>>>& chain_info)
{
    const int dof = num_layer * num_width;
    chain_num.assign(dof, 0);
    chain_info.assign(dof, std::vector<std::vector<int>>(max_chain_num, std::vector<int>(max_node_num, -1)));

    real_max_node_num = 0;
    for (int j = 0; j < num_layer; ++j)
    {
        for (int i = 0; i < num_width; ++i)
        {
            const int node_index = local_to_global_index(i, j, num_width);

            std::vector<int> neighbour_node_i;
            std::vector<int> neighbour_node_j;
            get_neighbour_node_list(i, j, neighbour_matrix[node_index], num_width, neighbour_node_i, neighbour_node_j);
            const int num_neighbour = neighbour_node_i.size();

            if (is_cross_linker[node_index] == 1)
            {

                chain_num[node_index] = num_neighbour;

                for (int k = 0; k < num_neighbour; ++k)
                {

                    std::vector<int> chain_node_list;
                    get_chain_node_list(i, j, neighbour_node_i[k], neighbour_node_j[k], num_width, entangle_orientation, is_cross_linker, chain_node_list);

                    chain_info[node_index][k] = chain_node_list;

                    real_max_node_num = std::max(real_max_node_num, int(chain_node_list.size()));
                }
            }

            else
            {

                chain_num[node_index] = 0;

                int neighbour_i = 0;
                int neighbour_j = 0;

                neighbour_i = i - 1;
                neighbour_j = j;
                std::vector<int> chain_node_list_left;
                get_chain_node_list(i, j, neighbour_i, neighbour_j, num_width, entangle_orientation, is_cross_linker, chain_node_list_left);

                neighbour_i = i;
                neighbour_j = j - 1;
                std::vector<int> chain_node_list_bottom;
                get_chain_node_list(i, j, neighbour_i, neighbour_j, num_width, entangle_orientation, is_cross_linker, chain_node_list_bottom);

                neighbour_i = i + 1;
                neighbour_j = j;
                std::vector<int> chain_node_list_right;
                get_chain_node_list(i, j, neighbour_i, neighbour_j, num_width, entangle_orientation, is_cross_linker, chain_node_list_right);

                neighbour_i = i;
                neighbour_j = j + 1;
                std::vector<int> chain_node_list_top;
                get_chain_node_list(i, j, neighbour_i, neighbour_j, num_width, entangle_orientation, is_cross_linker, chain_node_list_top);

                int node_l_index = local_to_global_index(i - 1, j, num_width);
                int node_r_index = local_to_global_index(i + 1, j, num_width);

                std::vector<int> neighbour_list = neighbour_matrix[node_index];

                auto it_l = std::find(neighbour_list.begin(), neighbour_list.end(), node_l_index);
                auto it_r = std::find(neighbour_list.begin(), neighbour_list.end(), node_r_index);

                if (entangle_orientation[node_index] == 0)
                {

                    {                                            

                        chain_node_list_left.erase(chain_node_list_left.begin());

                        std::reverse(chain_node_list_left.begin(), chain_node_list_left.end());

                        chain_node_list_left.insert(chain_node_list_left.end(), chain_node_list_bottom.begin(), chain_node_list_bottom.end());                
                    }

                    {

                        chain_node_list_right.erase(chain_node_list_right.begin());

                        std::reverse(chain_node_list_right.begin(), chain_node_list_right.end());

                        chain_node_list_right.insert(chain_node_list_right.end(), chain_node_list_top.begin(), chain_node_list_top.end());
                    }

                    if (it_l != neighbour_list.end() && it_r != neighbour_list.end())
                    {
                        chain_num[node_index] = 2;
                        chain_info[node_index][0] = chain_node_list_left;
                        chain_info[node_index][1] = chain_node_list_right;

                        real_max_node_num = std::max(real_max_node_num, int(chain_node_list_left.size()));
                        real_max_node_num = std::max(real_max_node_num, int(chain_node_list_right.size()));
                    }

                    else if (it_l == neighbour_list.end() && it_r != neighbour_list.end())
                    {
                        chain_num[node_index] = 1;
                        chain_info[node_index][0] = chain_node_list_right;

                        real_max_node_num = std::max(real_max_node_num, int(chain_node_list_right.size()));
                    }

                    else if (it_l != neighbour_list.end() && it_r == neighbour_list.end())
                    {
                        chain_num[node_index] = 1;
                        chain_info[node_index][0] = chain_node_list_left;

                        real_max_node_num = std::max(real_max_node_num, int(chain_node_list_left.size()));
                    }

                    else
                    {
                        chain_num[node_index] = 0;
                    }

                }

                else if (entangle_orientation[node_index] == 1)
                {

                    {

                        chain_node_list_left.erase(chain_node_list_left.begin());

                        std::reverse(chain_node_list_left.begin(), chain_node_list_left.end());

                        chain_node_list_left.insert(chain_node_list_left.end(), chain_node_list_top.begin(), chain_node_list_top.end());
                    }

                    {

                        chain_node_list_right.erase(chain_node_list_right.begin());

                        std::reverse(chain_node_list_right.begin(), chain_node_list_right.end());

                        chain_node_list_right.insert(chain_node_list_right.end(), chain_node_list_bottom.begin(), chain_node_list_bottom.end());
                    }

                    if (it_l != neighbour_list.end() && it_r != neighbour_list.end())
                    {
                        chain_num[node_index] = 2;
                        chain_info[node_index][0] = chain_node_list_left;
                        chain_info[node_index][1] = chain_node_list_right;

                        real_max_node_num = std::max(real_max_node_num, int(chain_node_list_left.size()));
                        real_max_node_num = std::max(real_max_node_num, int(chain_node_list_right.size()));
                    }

                    else if (it_l == neighbour_list.end() && it_r != neighbour_list.end())
                    {
                        chain_num[node_index] = 1;
                        chain_info[node_index][0] = chain_node_list_right;

                        real_max_node_num = std::max(real_max_node_num, int(chain_node_list_right.size()));
                    }

                    else if (it_l != neighbour_list.end() && it_r == neighbour_list.end())
                    {
                        chain_num[node_index] = 1;
                        chain_info[node_index][0] = chain_node_list_left;

                        real_max_node_num = std::max(real_max_node_num, int(chain_node_list_left.size()));
                    }

                    else
                    {
                        chain_num[node_index] = 0;
                    }
                }
            }
        }
    }

    if (real_max_node_num > max_node_num)
    {
        std::cerr << "Error: the chain node list is too long (with size " << real_max_node_num << ") in generate_entangle_network_entangle_orientation()" << std::endl;
        exit(1);
    }    

    for (int i = 0; i < dof; ++i)
    {
        for (int j = 0; j < chain_num[i]; ++j)
        {
            std::vector<int> node_in_chain = chain_info[i][j];
            node_in_chain = remove_const(node_in_chain, -1);

            chain_info[i][j] = node_in_chain;
        }
    }    
}

void compute_chain_init_len(
    int num_layer,
    int num_width,
    const std::vector<std::vector<double>>& x,
    const std::vector<std::vector<int>>& neighbour_matrix,
    const std::vector<int> & chain_num,
    const std::vector<std::vector<std::vector<int>>> & chain_info,
    std::vector<std::vector<double>> & dist_matrix_init,
    std::vector<std::vector<double>> & chain_init_len)
{   
    const int dof = num_layer * num_width;

    dist_matrix_init = compute_dist_matrix_simple(x, neighbour_matrix);
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
                        chain_len += dist_matrix_init[chain_info[i][j][k]][l];
                        break;
                    }
                }
            }

            chain_init_len[i][j] = chain_len;
        }
    }
}

void add_defect(
    int num_layer,
    int num_width,
    const std::vector<int> & is_cross_linker,
    const std::vector<int> & entangle_orientation,
    std::vector<std::vector<int>>& neighbour_matrix
    )
{

    int i = num_width / 2;
    int j = num_layer / 2;

    int neighbour_i = i;
    int neighbour_j = j - 1;

    int node_index = local_to_global_index(i, j, num_width);

    std::vector<int> chain_node_list;

    if (is_cross_linker[node_index] == 1)
    {        
        get_chain_node_list(i, j, neighbour_i, neighbour_j, num_width, entangle_orientation, is_cross_linker, chain_node_list);
    }
    else if (is_cross_linker[node_index] == 0)
    {

        neighbour_i = i - 1;
        neighbour_j = j;
        std::vector<int> chain_node_list_left;
        get_chain_node_list(i, j, neighbour_i, neighbour_j, num_width, entangle_orientation, is_cross_linker, chain_node_list_left);

        neighbour_i = i;
        neighbour_j = j - 1;
        std::vector<int> chain_node_list_bottom;
        get_chain_node_list(i, j, neighbour_i, neighbour_j, num_width, entangle_orientation, is_cross_linker, chain_node_list_bottom);

        neighbour_i = i + 1;
        neighbour_j = j;
        std::vector<int> chain_node_list_right;
        get_chain_node_list(i, j, neighbour_i, neighbour_j, num_width, entangle_orientation, is_cross_linker, chain_node_list_right);

        neighbour_i = i;
        neighbour_j = j + 1;
        std::vector<int> chain_node_list_top;
        get_chain_node_list(i, j, neighbour_i, neighbour_j, num_width, entangle_orientation, is_cross_linker, chain_node_list_top);

        if (entangle_orientation[node_index] == 0)
        {

            {                                            

                chain_node_list_left.erase(chain_node_list_left.begin());

                std::reverse(chain_node_list_left.begin(), chain_node_list_left.end());

                chain_node_list_left.insert(chain_node_list_left.end(), chain_node_list_bottom.begin(), chain_node_list_bottom.end());                
                chain_node_list = chain_node_list_left;
            }

        }

        else if (entangle_orientation[node_index] == 1)
        {

            {

                chain_node_list_right.erase(chain_node_list_right.begin());

                std::reverse(chain_node_list_right.begin(), chain_node_list_right.end());

                chain_node_list_right.insert(chain_node_list_right.end(), chain_node_list_bottom.begin(), chain_node_list_bottom.end());
                chain_node_list = chain_node_list_right;
            }
        }
    }

        int size_chain_node_list = chain_node_list.size();
        for (int k = 0; k < size_chain_node_list; ++k)
        {
            int node_index = chain_node_list[k];

            if (k == 0)
            {
                int remove_neighbour_node_index = chain_node_list[k+1];

                for (int l = 0; l < 4; ++l)
                {
                    if (neighbour_matrix[node_index][l] == remove_neighbour_node_index)
                    {
                        neighbour_matrix[node_index][l] = node_index;
                    }
                }
            }

            else if (k == size_chain_node_list - 1)
            {
                int remove_neighbour_node_index = chain_node_list[k-1];

                for (int l = 0; l < 4; ++l)
                {
                    if (neighbour_matrix[node_index][l] == remove_neighbour_node_index)
                    {
                        neighbour_matrix[node_index][l] = node_index;
                    }
                }
            }

            else
            {
                int remove_neighbour_node_index_1 = chain_node_list[k-1];
                int remove_neighbour_node_index_2 = chain_node_list[k+1];

                for (int l = 0; l < 4; ++l)
                {
                    if (neighbour_matrix[node_index][l] == remove_neighbour_node_index_1)
                    {
                        neighbour_matrix[node_index][l] = node_index;
                    }
                    if (neighbour_matrix[node_index][l] == remove_neighbour_node_index_2)
                    {
                        neighbour_matrix[node_index][l] = node_index;
                    }
                }
            }
        }
}

void generate_network_neighbour_matrix(
    int num_layer,
    int num_width,
    const std::vector<std::vector<double>>& x,
    std::vector<std::vector<int>>& neighbour_matrix,
    int is_fracture)
{
    assert((num_layer % 2 == 0) && (num_width % 2 == 0));
    const int half_num_layer = num_layer / 2;
    const int half_num_width = num_width / 2;

    const int dof = num_layer * num_width;

    for (int j = 0; j < num_layer; ++j)
    {
        for (int i = 0; i < num_width; ++i)
        {
            int node_index = local_to_global_index(i, j, num_width);

            int node_l_index = local_to_global_index(i - 1, j, num_width);      
            int node_r_index = local_to_global_index(i + 1, j, num_width);      
            int node_b_index = local_to_global_index(i, j - 1, num_width);      
            int node_t_index = local_to_global_index(i, j + 1, num_width);      

            std::vector<int> neighbour_list = { node_l_index, node_r_index, node_t_index, node_b_index };

            if (is_fracture != 0)
            {
                if (j == (half_num_layer - 1) && i <= (half_num_width - 1))
                {
                    neighbour_list = { node_l_index, node_r_index, node_index, node_b_index };
                }
                else if (j == (half_num_layer) && i <= (half_num_width - 1))
                {
                    neighbour_list = { node_l_index, node_r_index, node_t_index, node_index };
                }
            }

            for (int k = 0; k < 4; ++k)
            {
                int neighbour_index = neighbour_list[k];

                if (neighbour_index < 0 || neighbour_index >= dof)
                {
                    neighbour_index = node_index;
                }

                double d = dist_point(x[node_index], x[neighbour_index]);
                if (std::abs(d - 1) > 1e-10)
                {
                    neighbour_index = node_index;
                }
                neighbour_list[k] = neighbour_index;
            }

            neighbour_matrix[node_index] = neighbour_list;
        }
    }
}

void get_chain_node_list(int current_node_i, int current_node_j, 
        int neighbour_node_i, int neighbour_node_j,
        const int num_width,
        const std::vector<int> & entangle_orientation,
        const std::vector<int> & is_cross_linker,
        std::vector<int> & chain_node_list)
{
    chain_node_list.clear();

    int current_node_index = local_to_global_index(current_node_i, current_node_j, num_width);
    chain_node_list.push_back(current_node_index);

    int neighbour_node_index = local_to_global_index(neighbour_node_i, neighbour_node_j, num_width);
    chain_node_list.push_back(neighbour_node_index);

    int i = 0;
    while (true)
    {

        int last_node = chain_node_list.back();

        if (is_cross_linker[last_node] == 1)
        {
            return;
        }
        else
        {

            int next_neighbour_node_i = 0;
            int next_neighbour_node_j = 0;
            get_next_neighbour_node_coordinates(current_node_i, current_node_j, neighbour_node_i, neighbour_node_j, entangle_orientation[neighbour_node_index], next_neighbour_node_i, next_neighbour_node_j);

            int next_neighbour_node_index = local_to_global_index(next_neighbour_node_i, next_neighbour_node_j, num_width);
            chain_node_list.push_back(next_neighbour_node_index);

            current_node_i = neighbour_node_i;
            current_node_j = neighbour_node_j;

            neighbour_node_i = next_neighbour_node_i;
            neighbour_node_j = next_neighbour_node_j;
        }
        i++;

    }
}

void get_neighbour_node_list(const int i, const int j, const std::vector<int> & neighbour_list,
    const int num_width,
    std::vector<int> & neighbour_node_i, std::vector<int> & neighbour_node_j)
{    
    int node_l_index = local_to_global_index(i - 1, j, num_width);
    int node_r_index = local_to_global_index(i + 1, j, num_width);
    int node_b_index = local_to_global_index(i, j - 1, num_width);
    int node_t_index = local_to_global_index(i, j + 1, num_width);

    neighbour_node_i.clear();
    neighbour_node_j.clear();

    auto it_l = std::find(neighbour_list.begin(), neighbour_list.end(), node_l_index);
    if (it_l != neighbour_list.end())
    {
        neighbour_node_i.push_back(i - 1);
        neighbour_node_j.push_back(j);
    }

    auto it_r = std::find(neighbour_list.begin(), neighbour_list.end(), node_r_index);
    if (it_r != neighbour_list.end())
    {
        neighbour_node_i.push_back(i + 1);
        neighbour_node_j.push_back(j);
    }

    auto it_b = std::find(neighbour_list.begin(), neighbour_list.end(), node_b_index);
    if (it_b != neighbour_list.end())
    {
        neighbour_node_i.push_back(i);
        neighbour_node_j.push_back(j - 1);
    }

    auto it_t = std::find(neighbour_list.begin(), neighbour_list.end(), node_t_index);
    if (it_t != neighbour_list.end())
    {
        neighbour_node_i.push_back(i);
        neighbour_node_j.push_back(j + 1);
    }
}

void get_neighbour_node_orientation(const int current_node_i, const int current_node_j,
                                    const int neighbour_node_i, const int neighbour_node_j,
                                    std::string & neighbour_node_orientation)
{
    if (neighbour_node_i == current_node_i - 1)
    {
        neighbour_node_orientation = "l";
    }
    else if (neighbour_node_i == current_node_i + 1)
    {
        neighbour_node_orientation = "r";
    }
    else if (neighbour_node_j == current_node_j - 1)
    {
        neighbour_node_orientation = "b";
    }
    else if (neighbour_node_j == current_node_j + 1)
    {
        neighbour_node_orientation = "t";
    }
    else
    {
        std::cerr << "Error: invalid neighbour node orientation in get_neighbour_node_orientation()" << std::endl;
    }
}

void get_next_neighbour_node_coordinates(const int current_node_i, const int current_node_j,
                                         const int neighbour_node_i, const int neighbour_node_j,
                                         const int entangle_orientation_neighbour_node,
                                         int & next_neighbour_node_i, int & next_neighbour_node_j)
{
    assert(entangle_orientation_neighbour_node == 0 || entangle_orientation_neighbour_node == 1);

    std::string neighbour_node_orientation;
    get_neighbour_node_orientation(current_node_i, current_node_j, neighbour_node_i, neighbour_node_j, neighbour_node_orientation);

    if (neighbour_node_orientation == "l")
    {
        next_neighbour_node_i = current_node_i - 1;

        if (entangle_orientation_neighbour_node == 0)
        {            
            next_neighbour_node_j = current_node_j + 1;
        }

        else if (entangle_orientation_neighbour_node == 1)
        {
            next_neighbour_node_j = current_node_j - 1;
        }
    }

    else if (neighbour_node_orientation == "r")
    {
        next_neighbour_node_i = current_node_i + 1;

        if (entangle_orientation_neighbour_node == 0)
        {            
            next_neighbour_node_j = current_node_j - 1;
        }

        else if (entangle_orientation_neighbour_node == 1)
        {
            next_neighbour_node_j = current_node_j + 1;
        }
    }

    else if (neighbour_node_orientation == "b")
    {
        next_neighbour_node_j = current_node_j - 1;

        if (entangle_orientation_neighbour_node == 0)
        {
            next_neighbour_node_i = current_node_i + 1;
        }

        else if (entangle_orientation_neighbour_node == 1)
        {
            next_neighbour_node_i = current_node_i - 1;
        }
    }

    else if (neighbour_node_orientation == "t")
    {
        next_neighbour_node_j = current_node_j + 1;

        if (entangle_orientation_neighbour_node == 0)
        {
            next_neighbour_node_i = current_node_i - 1;
        }

        else if (entangle_orientation_neighbour_node == 1)
        {
            next_neighbour_node_i = current_node_i + 1;
        }
    }
}

void initial_cross_linker_random(const int num_layer, const int num_width,
                                const double prob_cross_linker,
                                const int random_seed,
                                const int is_fracture,
                                std::vector<int> & is_cross_linker)
{
    std::mt19937 gen(random_seed);
    std::uniform_real_distribution<> distr(0.0, 1.0);

    assert((num_layer % 2 == 0) && (num_width % 2 == 0));
    const int half_num_layer = num_layer / 2;
    const int half_num_width = num_width / 2;

    const int dof = num_layer * num_width;
    is_cross_linker.assign(dof, 0);

    for (int j = 0; j < num_layer; ++j)
    {
        for (int i = 0; i < num_width; ++i)
        {

            bool is_b = (j == 0);
            bool is_t = (j == (num_layer - 1));
            bool is_l = (i == 0);
            bool is_r = (i == (num_width - 1));

            int node_index = local_to_global_index(i, j, num_width);

            bool is_boundary = is_b || is_t || is_l || is_r;
            if (is_fracture != 0)
            {
                if (j == (half_num_layer - 1) && i <= (half_num_width - 1))
                {
                    is_boundary = true;
                }

                if (j == (half_num_layer) && i <= (half_num_width - 1))
                {
                    is_boundary = true;
                }
            }

            if (is_boundary)
            {
                is_cross_linker[node_index] = 1;
            }

            else
            {
                if (distr(gen) < prob_cross_linker)
                {
                    is_cross_linker[node_index] = 1;
                }
                else
                {
                    is_cross_linker[node_index] = 0;
                }
            }
        }
    }
}

void initial_cross_linker_random_min_gap(const int num_layer, const int num_width,
                                const double prob_cross_linker,
                                const int random_seed,
                                const int min_gap_x,
                                const int min_gap_y,
                                const int is_fracture,
                                std::vector<int> & is_cross_linker)
{
    const double scaled_prob_cross_linker = min_gap_x * min_gap_y * prob_cross_linker;

    if (scaled_prob_cross_linker > 1.0)
    {
        initial_cross_linker_random(num_layer, num_width, prob_cross_linker, random_seed, is_fracture, is_cross_linker);
        return;
    }

    assert((num_layer % 2 == 0) && (num_width % 2 == 0));
    const int half_num_layer = num_layer / 2;
    const int half_num_width = num_width / 2;

    std::mt19937 gen(random_seed);
    std::uniform_real_distribution<> distr(0.0, 1.0);

    std::uniform_int_distribution<int> distr_int_x(0, min_gap_x - 1);
    const int rand_int_x = distr_int_x(gen);
    std::uniform_int_distribution<int> distr_int_y(0, min_gap_y - 1);
    const int rand_int_y = distr_int_y(gen);

    const int dof = num_layer * num_width;
    is_cross_linker.assign(dof, 0);

    for (int j = 0; j < num_layer; ++j)
    {
        for (int i = 0; i < num_width; ++i)
        {

            bool is_b = (j == 0);
            bool is_t = (j == (num_layer - 1));
            bool is_l = (i == 0);
            bool is_r = (i == (num_width - 1));

            int node_index = local_to_global_index(i, j, num_width);

            bool is_boundary = is_b || is_t || is_l || is_r;
            if (is_fracture != 0)
            {
                if (j == (half_num_layer - 1) && i <= (half_num_width - 1))
                {
                    is_boundary = true;
                }

                if (j == (half_num_layer) && i <= (half_num_width - 1))
                {
                    is_boundary = true;
                }
            }

            bool is_crack_tip = false;
            if (is_fracture != 0)
            {
                if (j == (half_num_layer - 1) && i == (half_num_width))
                {
                    is_crack_tip = true;
                }

                if (j == (half_num_layer) && i == (half_num_width))
                {
                    is_crack_tip = true;
                }
            }

            if (is_boundary)
            {
                is_cross_linker[node_index] = 1;
            }

            else
            {
                if ((i % min_gap_x == rand_int_x) && (j % min_gap_y == rand_int_y))
                {
                    if (distr(gen) < scaled_prob_cross_linker)
                    {
                        is_cross_linker[node_index] = 1;
                    }
                }
            }

            if (is_crack_tip)
            {

            }            
        }
    }
}

void initial_cross_linker_random_min_gap_2(const int num_layer, const int num_width,
                                const double prob_cross_linker,
                                const int random_seed,
                                const int is_fracture,
                                std::vector<int> & is_cross_linker)
{
    const double scaled_prob_cross_linker = 2.0 * prob_cross_linker;

    if (prob_cross_linker >= 0.5)
    {
        initial_cross_linker_random(num_layer, num_width, prob_cross_linker, random_seed, is_fracture, is_cross_linker);
        return;
    }

    assert((num_layer % 2 == 0) && (num_width % 2 == 0));
    const int half_num_layer = num_layer / 2;
    const int half_num_width = num_width / 2;

    std::mt19937 gen(random_seed);
    std::uniform_real_distribution<> distr(0.0, 1.0);

    std::uniform_int_distribution<int> distr_int(0, 1);
    const int rand_int_0_1 = distr_int(gen);

    const int dof = num_layer * num_width;
    is_cross_linker.assign(dof, 0);

    for (int j = 0; j < num_layer; ++j)
    {
        for (int i = 0; i < num_width; ++i)
        {

            bool is_b = (j == 0);
            bool is_t = (j == (num_layer - 1));
            bool is_l = (i == 0);
            bool is_r = (i == (num_width - 1));

            int node_index = local_to_global_index(i, j, num_width);

            bool is_boundary = is_b || is_t || is_l || is_r;
            if (is_fracture != 0)
            {
                if (j == (half_num_layer - 1) && i <= (half_num_width - 1))
                {
                    is_boundary = true;
                }

                if (j == (half_num_layer) && i <= (half_num_width - 1))
                {
                    is_boundary = true;
                }
            }

            bool is_crack_tip = false;
            if (is_fracture != 0)
            {
                if (j == (half_num_layer - 1) && i == (half_num_width))
                {
                    is_crack_tip = true;
                }

                if (j == (half_num_layer) && i == (half_num_width))
                {
                    is_crack_tip = true;
                }
            }

            if (is_boundary)
            {
                is_cross_linker[node_index] = 1;
            }

            else
            {
                if (rand_int_0_1 == 0)
                {
                    if (((i % 2 == 1) && (j % 2 == 1)) || ((i % 2 == 0) && (j % 2 == 0)))
                    {
                        if (distr(gen) < scaled_prob_cross_linker)
                        {
                            is_cross_linker[node_index] = 1;
                        }
                    }
                }
                else
                {
                    if (((i % 2 == 1) && (j % 2 == 0)) || ((i % 2 == 0) && (j % 2 == 1)))
                    {
                        if (distr(gen) < scaled_prob_cross_linker)
                        {
                            is_cross_linker[node_index] = 1;
                        }
                    }

                }
            }

            if (is_crack_tip)
            {

            }            
        }
    }
}

void initial_cross_linker_quasi_random(const int num_layer, const int num_width,
                                const double prob_cross_linker,
                                const int random_seed,
                                std::vector<int> & is_cross_linker)
{
    const int dof = num_layer * num_width;
    is_cross_linker.assign(dof, 0);

    for (int j = 0; j < num_layer; ++j)
    {
        for (int i = 0; i < num_width; ++i)
        {

            bool is_b = (j == 0);
            bool is_t = (j == (num_layer - 1));
            bool is_l = (i == 0);
            bool is_r = (i == (num_width - 1));

            int node_index = local_to_global_index(i, j, num_width);

            if (is_b || is_t || is_l || is_r)
            {
                is_cross_linker[node_index] = 1;
            }
        }
    }

    const int x_lower_bound = 1;
    const int x_upper_bound = num_width - 2;
    const int y_lower_bound = 1;
    const int y_upper_bound = num_layer - 2;
    const int dof_inner = (num_layer - 2) * (num_width - 2);

    const int num_cross_linker_inner = std::ceil(prob_cross_linker * dof_inner);

    std::mt19937 gen(random_seed);
    std::uniform_real_distribution<> distr(0.0, 1.0);
    const double rand_num = distr(gen);

    std::vector<std::vector<int>> halton_seq_2d = generate_halton_seq_2d_integer(num_cross_linker_inner, x_lower_bound, x_upper_bound, y_lower_bound, y_upper_bound, rand_num);

    for (int i = 0; i < num_cross_linker_inner; ++i)
    {
        int x_index = halton_seq_2d[i][0];
        int y_index = halton_seq_2d[i][1];

        int node_index = local_to_global_index(x_index, y_index, num_width);
        is_cross_linker[node_index] = 1;
    }
}

double compute_cross_linker_ratio_inner_domain(int num_layer, int num_width, const std::vector<int> & is_cross_linker)
{
    int num_cross_linker = 0;
    int num_node = 0;
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
                num_node++;

                int node_index = local_to_global_index(i, j, num_width);
                if (is_cross_linker[node_index] == 1)
                {
                    num_cross_linker++;
                }
            }
        }
    }

    return num_cross_linker / (num_node + 0.0);
}

void initial_cross_linker_periodic_unit(const int num_layer, const int num_width,
                                        const int unit_size_x, const int unit_size_y,
                                        const int num_cross_linker_in_unit,
                                        const int is_fracture,
                                        std::vector<int> & is_cross_linker)
{
    assert((num_layer % 2 == 0) && (num_width % 2 == 0));
    const int half_num_layer = num_layer / 2;
    const int half_num_width = num_width / 2;

    for (int j = 0; j < num_layer; ++j)
    {
        for (int i = 0; i < num_width; ++i)
        {

            bool is_b = (j == 0);
            bool is_t = (j == (num_layer - 1));
            bool is_l = (i == 0);
            bool is_r = (i == (num_width - 1));

            bool is_boundary = is_b || is_t || is_l || is_r;
            if (is_fracture != 0)
            {
                if (j == (half_num_layer - 1) && i <= (half_num_width - 1))
                {
                    is_boundary = true;
                }

                if (j == (half_num_layer) && i <= (half_num_width - 1))
                {
                    is_boundary = true;
                }
            }

            bool is_crack_tip = false;
            if (is_fracture != 0)
            {
                if (j == (half_num_layer - 1) && i == (half_num_width))
                {
                    is_crack_tip = true;
                }

                if (j == (half_num_layer) && i == (half_num_width))
                {
                    is_crack_tip = true;
                }
            }

            if (!is_boundary)
            {

                if (unit_size_x == 2 && unit_size_y == 2)
                {

                    if (num_cross_linker_in_unit == 3)
                    {
                        if ((i % 2 == 1 && j % 2 == 1) 
                        || (i % 2 == 1 && j % 2 == 0)
                        || (i % 2 == 0 && j % 2 == 1))

                        {
                            is_cross_linker[i + j * num_width] = 1;
                        }
                        else
                        {
                            is_cross_linker[i + j * num_width] = 0;
                        }
                    }

                    else if (num_cross_linker_in_unit == 2)
                    {
                        if ((i % 2 == 0 && j % 2 == 1) 
                        || (i % 2 == 1 && j % 2 == 0))
                        {
                            is_cross_linker[i + j * num_width] = 1;
                        }
                        else
                        {
                            is_cross_linker[i + j * num_width] = 0;
                        }
                    }

                    else if (num_cross_linker_in_unit == 1)
                    {
                        if (i % 2 == 0 && j % 2 == 1)

                        {
                            is_cross_linker[i + j * num_width] = 1;
                        }
                        else
                        {
                            is_cross_linker[i + j * num_width] = 0;
                        }
                    }

                    else if (num_cross_linker_in_unit == 0)
                    {
                        is_cross_linker[i + j * num_width] = 0;
                    }
                }
            }
            else
            {
                is_cross_linker[i + j * num_width] = 1;
            }

            if (is_crack_tip)
            {
                is_cross_linker[i + j * num_width] = 0;
            }
        }
    }    
}