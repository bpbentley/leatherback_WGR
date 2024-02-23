#!/bin/bash
#SBATCH -J mosdepth
#SBATCH -o ./logs/mosdepth.log
#SBATCH -e ./logs/mosdepth.err
#SBATCH -p cpu-long
#SBATCH -t 12:00:00
#SBATCH -c 12  # Number of Cores per Task
#SBATCH --mem=20000  # Requested Memory

module load miniconda/4.11.0

conda activate turtles

mosdepth -n --by 10000 ./data/alignments/dc_33129_DR_RG_IR_SORT.bam

conda deactivate