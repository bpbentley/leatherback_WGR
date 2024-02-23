#!/bin/bash
#SBATCH -J PCAngsd_WGR_DerCor_filt
#SBATCH -o ./logs/PCAngsd/PCAngsd_WGR_DerCor_filt_NGSAdmix_%a.log
#SBATCH -e ./logs/PCAngsd/PCAngsd_WGR_DerCor_filt_NGSAdmix_%a.err
#SBATCH -p cpu-long
#SBATCH -t 12:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=64000  # Requested Memory
#SBATCH --array=2-15

module load miniconda/22.11.1-1
module load angsd/0.935

INDIR=/nese/meclab/Blair/WGR_DerCor/ANGSD/for_PCAngsd
OUTDIR=/nese/meclab/Blair/WGR_DerCor/PCAngsd
NGSDIR=/nese/meclab/Shared/bin/NGSadmix
k=${SLURM_ARRAY_TASK_ID}

conda activate pcangsd

pcangsd -b $INDIR/DerCor_filt3.beagle.gz --admix -t 16 -e $k -o $OUTDIR/DerCor_filt_Admix_$k

conda deactivate


OUTDIR=/nese/meclab/Blair/WGR_DerCor/NGSadmix

${NGSDIR}/NGSadmix -likes ${INDIR}/DerCor_filt3.beagle.gz -K ${SLURM_ARRAY_TASK_ID} -minMaf 0.05 -minInd 20 -maxiter 5000 -P 16 -outfiles ${OUTDIR}/DerCor_filt_${k}