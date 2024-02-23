#!/bin/bash
#SBATCH -J wpac_PCAngsd_data_generation_beagle
#SBATCH -o ./logs/ANGSD_PCangsd/wpac_PCAngsd_data_generation_filt_beagle_%a_%A.log
#SBATCH -e ./logs/ANGSD_PCangsd/wpac_PCAngsd_data_generation_filt_beagle_%a_%A.err
#SBATCH -p cpu-long
#SBATCH -t 72:00:00
#SBATCH -c 20  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=96000  # Requested Memory
#SBATCH --array=1-28

OUTDIR=/nese/meclab/Blair/WGR_DerCor/ANGSD/for_PCAngsd/wpac
CHR=SUPER_${SLURM_ARRAY_TASK_ID}
REFDIR=/nese/meclab/Shared/reference_genomes/ST_reference_genomes/rDerCor1_202105
REF=rDerCor1.pri.cur.20210524.fasta

module load angsd/0.935

angsd -GL 1 -out $OUTDIR/DerCor_wpac_${CHR} -nThreads 20 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -bam ./wpac_bams.txt -r $CHR -ref $REFDIR/$REF