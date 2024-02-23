#!/bin/bash
#SBATCH -J downsample_WGR_2x
#SBATCH -o ./logs/downsample/downsample_WGR_2x_%a_%A.log
#SBATCH -e ./logs/downsample/downsample_WGR_2x_%a_%A.err
#SBATCH -p cpu-long
#SBATCH -t 24:00:00
#SBATCH -c 12  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=32000  # Requested Memory
#SBATCH --array=2-110

module load picard/3.1.0

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p downsample_list_2X.txt | cut -f1)
DS=$(sed -n ${SLURM_ARRAY_TASK_ID}p downsample_list_2X.txt | cut -f3)
INDIR=/nese/meclab/Blair/WGR_DerCor/data/alignments
OUTDIR=/nese/meclab/Blair/WGR_DerCor/data/alignments/DS_2X
INBAM=${SAMPLE}_DR_RG_IR_SORT.bam
OUTBAM=${SAMPLE}_DR_RG_IR_SORT_DS2X.bam

picard DownsampleSam -I ${INDIR}/${INBAM} -O ${OUTDIR}/${OUTBAM} --CREATE_INDEX true -P ${DS}