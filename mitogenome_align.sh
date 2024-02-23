#!/bin/bash
#SBATCH -J mtDNA_align
#SBATCH -o ./logs/mito_align/mito_align_%A_%a.log
#SBATCH -e ./logs/mito_align/mito_align_%A_%a.err
#SBATCH -p cpu-long
#SBATCH -t 48:00:00
#SBATCH -c 12  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=16000  # Requested Memory
#SBATCH --array=5-11,15-19,46-48,50-62

module load bowtie/2.4.5
module load samtools/1.14+py3.8.12

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./all_IDs.txt)
CONCAT=/nese/meclab/Shared/raw_reads/DerCor_WGR_concat/concat
REFDIR=/nese/meclab/Shared/reference_genomes/ST_reference_genomes/rDerCor1_202105
REF=rDerCor1.MT.20190820
OUTDIR=/nese/meclab/Blair/WGR_DerCor/mtDNA
METDIR=/nese/meclab/Blair/WGR_DerCor/metrics/alignment_metrics/mtDNA

#bowtie2 -p 16 -1 $CONCAT/${SAMPLE}_1.fq.gz -2 $CONCAT/${SAMPLE}_2.fq.gz \
# -x $REFDIR/$REF -S ${OUTDIR}/${SAMPLE}_mitogenome.sam 2> ${METDIR}/${SAMPLE}_mitogenome_alignment_stats.txt | 
 
samtools view -@ 16 -S -b ${OUTDIR}/sams/${SAMPLE}_mitogenome.sam > ${OUTDIR}/${SAMPLE}_mitogenome.bam
rm ${OUTDIR}/sams/${SAMPLE}_mitogenome.sam