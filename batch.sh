#!/bin/bash
#SBATCH -N 1
#SBATCH -p testing
#SBATCH -t 4:00:00
#SBATCH --ntasks-per-node 1
 
g++ main_filter.cpp -o a.out -std=gnu++11
./a.out

 

wait # IMPORTANT: wait for all to finish or get kill
