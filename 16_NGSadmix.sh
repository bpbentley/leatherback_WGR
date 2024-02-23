#!/bin/bash
#SBATCH -J NGSadmix_DerCor_filt
#SBATCH -o ./logs/NGSadmix/NGSadmix_DerCor_filt_%a.log
#SBATCH -e ./logs/NGSadmix/NGSadmix_DerCor_filt_%a.err
#SBATCH -p cpu
#SBATCH -t 04:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=96000  # Requested Memory
#SBATCH --array=2-15

NGSDIR=/nese/meclab/Shared/bin/NGSadmix
INDIR=
OUTDIR=

${NGSDIR}/NGSadmix -likes ${INDIR}/DerCor_filt.beagle.gz -K ${SLURM_ARRAY_TASK_ID} -minMaf 0.05 -minInd 20 -outfiles ${OUTDIR}/DerCor_filt