#!/bin/bash
#SBATCH -N 1
#SBATCH -p defq
#SBATCH -t 10:00:00
#SBATCH --ntasks-per-node 1
 
g++ main-v3.0-grid.cpp -o grid.out -std=gnu++11
./grid.out

 

wait # IMPORTANT: wait for all to finish or get kill
