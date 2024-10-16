# Kidney Exchange Optimization - Group 7

This repository contains code and instance files related to kidney exchange optimization problems, with a focus on both exact and heuristic methods. It implements different approaches for solving kidney exchange problems, including the branch-and-bound algorithm and cycle-chain deactivation strategies.

## Table of Contents

1. [Project Structure](#project-structure)
2. [Getting Started](#getting-started)
3. [Exact Methods](#exact-methods)
4. [Heuristic Methods](#heuristic-methods)
5. [Instance Files](#instance-files)
6. [Usage](#usage)
7. [License](#license)

## Project Structure

The repository is structured as follows:

- **Exact Methods/**: Contains implementations of exact optimization methods.
  - `branch_and_bound/`: Branch-and-bound algorithm implemented in Jupyter notebooks.
    - `branch_and_bound_BFS.ipynb`: Implementation of branch-and-bound using breadth-first search approach for finding cycles and chains.
    - `branch_and_bound_recursive.ipynb`: Implementation of branch-and-bound using recursive approach for finding cycles and chains.
  - `cycle_chain_deactivation/`: Contains Python scripts for solving the KEP using the cycle-chain deactivation method.
    - `allocation.py`: Standard allocation file containing a variety of classes.
    - `allocation_generalized.py`: Generalized allocation file containing variety of classes and functions, including the BFS algorithm for finding cycles and chains.
    - `cycle_chain_deactivation.py`: Original cycle and chain deactivation algorithm wiyh fixed cycle and chain lengths.
    - `cycle_chain_deactivation_generalized.py`: Generalized version of the cycle-chain deactivation algorithm.
    - `run.ipynb`: Jupyter notebook for running cycle-chain deactivation on one file.
    - `run.py`: Python script to run kidney exchange optimization using normal or generalized methods. It accepts input/output files and options for cycle and chain lengths.
    - `run.sh`: Shell script that automates running instances with options to skip files, use generalized method, and set cycle/chain lengths.

- **Heuristic Methods/**: Contains heuristic algorithms for solving kidney exchange problems.
  - `allocation.py`: Allocation file for heuristics (removed Gurobi aspects).
  - `notebook.ipynb`: Jupyter notebook for testing and running heuristic method.
  - `run.py`: Python script to run the full heuristic approach for multiple files.

- **Instance Files/**: Contains input data files for the optimization methods.
  - Various files representing different kidney exchange instances, with and without non-directed donors (NDDs).

## Getting Started

To get started, clone the repository and set up the necessary dependencies. You will need Python installed, along with the required libraries specified in the project.

### Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/KidneyExchangeOptimization.git
cd KidneyExchangeOptimization
