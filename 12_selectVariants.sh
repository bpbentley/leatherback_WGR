#!/bin/bash
#SBATCH -J WGR_SelectVariants
#SBATCH -o ./logs/GenomicsBD/SelectVariants_%a.log
#SBATCH -e ./logs/GenomicsBD/SelectVariants_%a.err
#SBATCH -p gpu
#SBATCH -t 04:00:00
#SBATCH -c 12  # Number of Cores per Task
#SBATCH --mem=64000  # Requested Memory
#SBATCH --array=2-28

module load gatk/4.2.2.0+py3.8.12
module load htslib/1.12
module load samtools/1.14+py3.8.12

REFDIR=/nese/meclab/Shared/reference_genomes/ST_reference_genomes/rDerCor1_202105
REFERENCE=rDerCor1.pri.cur.20210524.fasta
INDIR=/nese/meclab/Blair/WGR_DerCor/snpEff/vcf
OUTDIR=/nese/meclab/Blair/WGR_DerCor/snpEff/tvcf

# Adjust NAME to whatever prefix you want your output files to have
VCF=WGR_DerCor_SUPER_${SLURM_ARRAY_TASK_ID}.vcf # unique name needed if running multiple pipelines on different data, as files will be written to same directory, and will over-write if not unique.
TVCF=DerCor_SV_SUPER_${SLURM_ARRAY_TASK_ID}.vcf
CHR=SUPER_${SLURM_ARRAY_TASK_ID}


gatk \
    SelectVariants \
    --java-options "-Xms60G -Xmx90G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -R ${REFDIR}/${REFERENCE} \
    --select-type-to-include SNP \
    --select-type-to-include INDEL \
    -L ${CHR} \
    -V ${INDIR}/${VCF} \
    -O ${OUTDIR}/${TVCF}