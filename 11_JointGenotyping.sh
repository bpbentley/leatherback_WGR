#!/bin/bash
#SBATCH -J WGR_JointGenotype
#SBATCH -o ./logs/GenomicsBD/WGR_JointGenotype_%a.log
#SBATCH -e ./logs/GenomicsBD/WGR_JointGenotype_%a.err
#SBATCH -p cpu-long
#SBATCH -t 72:00:00
#SBATCH -c 12  # Number of Cores per Task
#SBATCH --mem=64000  # Requested Memory
#SBATCH --array=2-28

module load gatk/4.2.2.0+py3.8.12
module load htslib/1.12
module load samtools/1.14+py3.8.12

REFDIR=/nese/meclab/Shared/reference_genomes/ST_reference_genomes/rDerCor1_202105
REF=rDerCor1.pri.cur.20210524.fasta
OUTDIR=/nese/meclab/Blair/WGR_DerCor/snpEff/vcf


gatk --java-options "-Xmx40g -Xms40g" GenotypeGVCFs \
      -R ${REFDIR}/${REF} \
      -V gendb:///nese/meclab/Blair/WGR_DerCor/snpEff/SUPER_${SLURM_ARRAY_TASK_ID} \
      -O ${OUTDIR}/WGR_DerCor_SUPER_${SLURM_ARRAY_TASK_ID}.vcf