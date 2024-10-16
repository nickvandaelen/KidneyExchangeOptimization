# Kidney Exchange Optimization - Group 7

This repository contains code and instance files related to kidney exchange optimization problems. It implements different approaches for solving kidney exchange problems, including the branch-and-bound algorithm, cycle-chain deactivation strategies and a variety of heuristic methods.

## Getting Started

To get started, clone the repository and set up the necessary dependencies. You will need Python installed, along with the required libraries specified in the project.

### Installation

Clone the repository:

```bash 
git clone https://github.com/nickvandaelen/KidneyExchangeOptimization.git
cd KidneyExchangeOptimization
```

## Project Structure

The repository is structured as follows:

- **Exact Methods/**: Contains implementations of exact optimization methods.
  - `branch_and_bound/`: Branch-and-bound algorithm implemented in Jupyter notebooks.
    - `branch_and_bound_BFS.ipynb`: Implementation of branch-and-bound using breadth-first search approach for finding cycles and chains.
    - `branch_and_bound_recursive.ipynb`: Implementation of branch-and-bound using recursive approach for finding cycles and chains.
  - `cycle_chain_deactivation/`: Contains Python scripts for solving the KEP using the cycle-chain deactivation method.
    - `allocation.py`: Standard allocation file containing a variety of classes.
    - `allocation_generalized.py`: Generalized allocation file containing variety of classes and functions, including the BFS algorithm for finding cycles and chains.
    - `cycle_chain_deactivation.py`: Original cycle and chain deactivation algorithm with fixed cycle and chain lengths.
    - `cycle_chain_deactivation_generalized.py`: Generalized version of the cycle-chain deactivation algorithm.
    - `run.ipynb`: Jupyter notebook for running cycle-chain deactivation on one file.
    - `run.py`: Python script to run kidney exchange optimization using normal or generalized methods. It accepts input/output files and options for cycle and chain lengths.
    - `run.sh`: Shell script that automates running instances with options to skip files, use generalized method, and set cycle/chain lengths.
      #### Run a specific instance:
      ```bash
       ./run.sh -f Delorme_1000_NDD_Unit_0.txt
      ```
      
      #### Run all instances:
      ```bash
      ./run.sh -a
      ```
      
      #### Run with generalized version:
      ```bash
      ./run.sh -f Delorme_1000_NDD_Unit_0.txt -g -c 4 -h 4
      ```
      ```bash
      ./run.sh -a -g -c 4 -h 4
      ```
      
      #### Run with skip patterns:
      ```bash
      ./run.sh -f Delorme_1000_NDD_Unit_0.txt -g -c 4 -h 4 -s "1000" # skips all files with 1000 nodes
      ```
      ```bash
      ./run.sh -a -g -c 4 -h 4 -s "Saidman" # skips all Saidman files
      ```
      ```bash
      ./run.sh -a -g -c 4 -h 4 -s "Delorme_1000" -s "Saidman_500" # skips all Delorme files with 1000 nodes and Saidman files with 500 nodes
      ```

- **Heuristic Methods/**: Contains heuristic algorithms for solving kidney exchange problems.
  - `allocation.py`: Allocation file for heuristics (removed Gurobi aspects).
  - `notebook.ipynb`: Jupyter notebook for testing and running heuristic method.
  - `run.py`: Python script to run the full heuristic approach for multiple files.

- **Instance Files/**: Contains input data files for the optimization methods.
  - Various files representing different kidney exchange instances.
