#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=00:10:00
#SBATCH --partition=cpu_short
#SBATCH --job-name=regression_fast
#SBATCH --output=%x.o%j
#SBATCH --export=PATH,LD_LIBRARY_PATH,FONTCONFIG_PATH,HOME,USER

module purge        

module load openmpi/4.1.5/gcc-13.2.0
module load anaconda3/2024.06/gcc-13.2.0

source activate my_r_env

cd basic_model

srun  Rscript toy.r