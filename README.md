# Entangled Network Simulation

This repository contains C++ and MATLAB codes for simulating and visualizing the mechanics of entangled networks.
For more details on the algorithm and implementation, please refer to our paper:

[J. Huang, J. Liu, and S. Lin. Topological Mechanics of Entangled Networks. *arXiv preprint arXiv:2509.17813*. (2025)](https://arxiv.org/abs/2509.17813)

---

## Environment

- **Compiler:** GCC 11
- **Platform:** Linux / macOS  
- **Dependencies:** Standard C++ libraries and MATLAB (for plotting)

---

## Compile and Run

### 1. Compile

```bash
cd code_cpp
make
````

### 2. Run Simulation

```bash
./main 20 80 10 5.0 40 1e-2 100000 1e-1 1e-1 5 0.0 2 2 2 0 0.0 1e-3 0 0.5 100 1 1234
```
This command runs the simulation with the specified parameters. The computation typically takes a few minutes to complete, depending on the hardware configuration.

After the run finishes, the program will automatically generate:
- Log files (*.txt) containing runtime information and diagnostic data.
- Binary data files (*.bin) storing numerical results.

You can then proceed to the next step to visualize and analyze the results in MATLAB.

### Example Input Parameters

| Argument                       | Description                             | Example Value |
| ------------------------------ | --------------------------------------- | ------------- |
| num_layer                      | Number of layers                        | 20            |
| num_width                      | Number of nodes per layer               | 80            |
| max_node_num                   | Maximum number of nodes                 | 10            |
| total_stretch                  | Total applied stretch                   | 5.0           |
| num_sampling_stretch           | Number of sampled stretch points        | 40            |
| error_tol                      | Error tolerance                         | 1e-2          |
| max_time_step                  | Maximum number of time steps            | 100000        |
| dt_no_entangle                 | Time step (no entanglement)             | 1e-1          |
| dt_entangle                    | Time step (entangled)                   | 1e-1          |
| num_decay_dt                   | Number of decay steps                   | 5             |
| rand_loc                       | Random location flag                    | 0.0           |
| unit_size_x                    | Unit size in x                          | 2             |
| unit_size_y                    | Unit size in y                          | 2             |
| num_cross_linker_in_unit       | Number of crosslinkers per unit         | 2             |
| is_random_cross_linker         | Use random crosslinkers (0/1)           | 0             |
| prob_cross_linker              | Probability of crosslinker              | 0.0           |
| error_tol_energy_stretch       | Energy–stretch tolerance                | 1e-3          |
| is_random_entangle_orientation | Random entanglement orientation (0/1)   | 0             |
| prob_entangle_orientation      | Probability of entanglement orientation | 0.5           |
| print_time_interval            | Print interval (steps)                  | 100           |
| is_fracture                    | Enable fracture (0/1)                   | 1             |
| rand_seed                      | Random seed                             | 1234          |

---

## Visualization

To generate plots from the simulation results:

```bash
cd code_matlab_plot
run main.m
```

This script will automatically produce several figures visualizing the networks in the undeformed and deformed states.

---

## Repository Structure

```
.
├── code_cpp/             # C++ source files and Makefile
│   ├── main.cpp
│   ├── ...
│   └── Makefile
├── code_matlab_plot/     # MATLAB plotting scripts
│   ├── main.m
│   └── ...
└── README.md             # This file
```

---

**Contact:**
For questions or feedback, please reach out to the authors via GitHub Issues or email (huangjt@udel.edu).

```
