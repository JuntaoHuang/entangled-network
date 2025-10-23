
clc; clear; close all;

path = '../code_cpp/';

num_layer = 20; num_width = 80; maxnode = 10;
max_stretch = 5; sampling = 40;
maxstp = 100000;
rand_percent = 0.0;
unit_size_x = 2; unit_size_y = 2; num_cross_linker_in_unit = 2;
is_random_cross_linker = 0;
prob_cross_linker = 0.0;
print_step = 100;
is_fracture = 1;
seed = 1234;

filename = 'lay_' + string(num_layer) ...
        + '_wid_' + string(num_width) ...
        + '_maxnode_' + string(maxnode) ...
        + '_stret_' + string(max_stretch) ...
        + '_samp_' + string(sampling) ...
        + '_tol_2_maxstp_' + string(maxstp) ...
        + '_dtne_1_dte_1_ntdcy_5_' ...
        + 'randloc_' + string(rand_percent * 100) ...
        + '_unitx_' + string(unit_size_x) ...
        + '_unity_' + string(unit_size_y) ...
        + '_ncross_' + string(num_cross_linker_in_unit) ...
        + '_randcross_' + string(is_random_cross_linker) ...
        + '_pcross_' + string(prob_cross_linker * 100) ...
        + '_energytol_3_randent_0_pent_50' ...
        + '_print_' + string(print_step) ...
        + '_frac_' + string(is_fracture) ...
        + '_seed_' + string(seed) ...
        + '_cpp';
filename = char(filename);

x_init_cpp = loadVector2DToArrayDouble([path 'x_init_' filename '.bin']);

x_entangle_cpp = loadVector2DToArrayDouble([path 'x_entangle_crack_' filename '.bin']);
x_no_entangle_cpp = loadVector2DToArrayDouble([path 'x_no_entangle_crack_' filename '.bin']);

chain_num_cpp = loadVector1DToArrayInt([path 'chain_num_' filename '.bin']);
chain_init_len_cpp = loadVector2DToArrayDouble([path 'chain_init_len_' filename '.bin']);
chain_info_cpp = loadVector3DToArrayInt([path 'chain_info_' filename '.bin']);
is_cross_linker_cpp = loadVector1DToArrayInt([path 'is_cross_linker_' filename '.bin']);
neighbour_matrix_cpp = loadVector2DToArrayInt([path 'neighbour_matrix_' filename '.bin']);

%% plot network
% initial state
fig_name = ['network_init_no_entangle_' filename];
is_entangle = false;
plot_network(x_init_cpp, neighbour_matrix_cpp, is_cross_linker_cpp, is_entangle, fig_name);

fig_name = ['network_init_entangle_' filename];
is_entangle = true;
plot_network(x_init_cpp, neighbour_matrix_cpp, is_cross_linker_cpp, is_entangle, fig_name);

% after deformation
fig_name = ['network_no_entangle_stretch_' filename];
plot_network_stretch_no_entangle(x_no_entangle_cpp, neighbour_matrix_cpp, x_init_cpp, fig_name);

fig_name = ['network_entangle_stretch_' filename];
is_entangle = true;
plot_network_stretch_entangle(x_entangle_cpp, neighbour_matrix_cpp, chain_num_cpp, chain_init_len_cpp, chain_info_cpp, is_cross_linker_cpp, fig_name);