#!/bin/bash
#SBATCH -N 1
#SBATCH -p defq
#SBATCH -t 5:00:00
#SBATCH --ntasks-per-node 1
 
g++ main_nongrid.cpp -o b.out -std=gnu++11
./b.out

 

wait # IMPORTANT: wait for all to finish or get kill
