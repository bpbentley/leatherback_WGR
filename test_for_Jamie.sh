#!/bin/bash -l

#SBATCH --time 15:00:00
#SBATCH -p cpu-long
#SBATCH -J "angsdVCF"
#SBATCH -o ./logs/log_%j%x
#SBATCH -c 6
#SBATCH --mem=96000


module load angsd/0.935

nInd=96

minInd=1

angsd -bam ./filtered_alignments.txt -out restrict_2019HIRAD -P 5 -minQ 20 -minMapQ 10 -minInd $minInd -dobcf 1 -doPost 1 -postCutoff 0.95 -doGeno 3 -GL 1 -doMajorMinor 1 -doMaf 2
