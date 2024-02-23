#!/bin/bash
#SBATCH -J WGR_GenomicsBD
#SBATCH -o ./logs/GenomicsBD/WGR_GenomicsBD_%a.log
#SBATCH -e ./logs/GenomicsBD/WGR_GenomicsBD_%a.err
#SBATCH -p cpu-long
#SBATCH -t 72:00:00
#SBATCH -c 12  # Number of Cores per Task
#SBATCH --mem=20000  # Requested Memory
#SBATCH --array=2-28 #(N=40)

module load gatk/4.2.2.0+py3.8.12
module load htslib/1.12
module load samtools/1.14+py3.8.12

OUTDIR=/project/uma_lisa_komoroske/Blair/WGR_for_genome/snpEff/CheMyd/exons/

gatk --java-options "-Xmx40g -Xms40g" GenomicsDBImport \
      --sample-name-map cohort.sample_map \
      --genomicsdb-workspace-path /nese/meclab/Blair/WGR_DerCor/snpEff/SUPER_${SLURM_ARRAY_TASK_ID} \
      -L SUPER_${SLURM_ARRAY_TASK_ID}