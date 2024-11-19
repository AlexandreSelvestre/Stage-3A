#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH --time=4:00:00
#SBATCH --partition=cpu_med
#SBATCH --job-name=regression_fast
#SBATCH --output=%x.o%j
#SBATCH --export=PATH,LD_LIBRARY_PATH,FONTCONFIG_PATH,HOME,USER

module purge        
module load libxml2/2.9.9/gcc-9.2.0     
module load zlib/1.2.11/gcc-11.2.0     
module load openssl/1.1.1l/gcc-11.2.0   
module load gsl/2.6/gcc-9.2.0
module load zlib/1.2.9/gcc-9.2.0
module load libpng/1.6.37/gcc-9.2.0
module load zstd/1.4.3/gcc-9.2.0
module load openmpi/4.1.5/gcc-13.2.0
module load gcc/13.2.0/gcc-4.8.5
module load anaconda3/2024.06/gcc-13.2.0

source activate my_r_env

cd ../basic_model

srun  Rscript run.r
