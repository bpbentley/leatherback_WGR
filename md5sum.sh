#!/bin/bash
#SBATCH -J DerCor_WGR_md5sum # Job name
#SBATCH -c 4  # Number of Cores per Task
#SBATCH --mem=8192  # Requested Memory
#SBATCH -p cpu  # Partition
#SBATCH -t 10:00:00  # Job time limit
#SBATCH -o ./logs/md5sums/md5sum.%j.out  # %j = job ID
#SBATCH -e ./logs/md5sums/md5sum.%j.err  # %j = job ID
#SBATCH --array=1-167%50

FILE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./all_alignments.txt)
SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./all_samples.txt)

md5sum ${FILE} > md5sums/${SAMPLE}_md5.txt