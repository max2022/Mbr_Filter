#!/bin/bash
#SBATCH -N 1
#SBATCH -p defq
#SBATCH -t 40:00:00
#SBATCH --ntasks-per-node 1
 
g++ main_test_1.cpp -o c.out -std=gnu++11
./c.out

 

wait # IMPORTANT: wait for all to finish or get kill
