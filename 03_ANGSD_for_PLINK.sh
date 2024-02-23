#!/bin/bash
#SBATCH -J 2X_PLINK_data_generation
#SBATCH -o ./logs/ANGSD_PLINK/2X_PLINK_data_generation_filt_%a_%A.log
#SBATCH -e ./logs/ANGSD_PLINK/2X_PLINK_data_generation_filt_%a_%A.err
#SBATCH -p cpu-long
#SBATCH -t 72:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=96000  # Requested Memory
#SBATCH --array=1-28


module load samtools/1.14
module load angsd/0.935
module load bedtools2/2.30.0+py3.8.12

REFDIR=/nese/meclab/Shared/reference_genomes/ST_reference_genomes/rDerCor1_202105
REF=rDerCor1.pri.cur.20210524.fasta
NUM=$(printf ${SLURM_ARRAY_TASK_ID})
CHR=SUPER_${NUM}
OUTDIR=/nese/meclab/Blair/WGR_DerCor/ANGSD/for_PLINK

angsd -bam ./ds_bams.txt -out ${OUTDIR}/DerCor_2X_${CHR} -doPlink 2 -doGeno -4 -doPost 1 -doMajorMinor 1 -GL 1 -doCounts 1 -doMaf 2 -postCutoff 0.95 -SNP_pval 1e-6 -uniqueonly 1 -remove_bads 1 -C 50 -baq 1 -r $CHR -ref $REFDIR/$REF -nThreads 20