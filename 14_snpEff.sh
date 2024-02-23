#!/bin/bash
#SBATCH -J WGR_snpEff
#SBATCH -o ./logs/snpEff/WGR_snpEff_%a_%A.log
#SBATCH -e ./logs/snpEff/WGR_snpEff_%a_%A.err
#SBATCH -p cpu
#SBATCH -t 01:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=32000  # Requested Memory
#SBATCH --array=2-40

module load bcftools/1.15
module load snpeff/2017-11-24

config=/nese/meclab/Blair/WGR_DerCor/snpEff/rDerCor1.pri.cur.20210524/snpEff.config
OUTDIR=/nese/meclab/Blair/WGR_DerCor/snpEff/snpEff_out
INDIR=/nese/meclab/Blair/WGR_DerCor/snpEff/fvcf
SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./high_coverage_IDs.txt)

bcftools view -c1 -Ob -s ${SAMPLE} -o $INDIR/${SAMPLE}.vcf ${INDIR}/DerCor_all_recode.vcf

snpEff eff rDerCor1.pri.cur.20210524 -nodownload -c ${config} -s $OUTDIR/${SAMPLE} $INDIR/${SAMPLE}.vcf > $INDIR/${SAMPLE}.snpEff.vcf