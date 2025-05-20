#!/bin/bash -l
#SBATCH --account jriou_ofsp_surveillance 
#SBATCH --mail-user julien.riou@unisante.ch
#SBATCH --job-name pfas2
#SBATCH --partition cpu
#SBATCH --cpus-per-task 4
#SBATCH --mem 10G 
#SBATCH --time 96:00:00 

module load r-light 

Rscript index.R 

