#!/bin/bash
#BSUB -J "DerCor_WGR_align_process[1-167]"
#BSUB -o ./logs/02_align_process/DerCor_WGR_generateBAMs_%I.log
#BSUB -e ./logs/02_align_process/DerCor_WGR_generateBAMs_%I.err
#BSUB -q long
#BSUB -W 120:00
#BSUB -R rusage[mem=6000]
#BSUB -R span[hosts=1]
#BSUB -n 16

module load bwa/0.7.17
module load samtools/1.9
module load python3/3.8.2
module load fastqc/0.11.5
module load trimmomatic/0.39
module load picard/2.23.3
module load htslib/1.9
module load anaconda3/2019.03
module load gatk/3.5

REFDIR=/project/uma_lisa_komoroske/Blair/refs/rDerCor1_20210524
REF=rDerCor1.pri.cur.20210524.fasta
FQDIR=/project/uma_lisa_komoroske/Blair/DerCor_WGR_Final/concat_files
SAMPLE=$(sed -n ${LSB_JOBINDEX}'p' ./all_samples.txt)
OUTDIR=/project/uma_lisa_komoroske/Blair/DerCor_WGR_Final/alignments
METDIR=/project/uma_lisa_komoroske/Blair/DerCor_WGR_Final/metrics
TRIMDIR=/project/uma_lisa_komoroske/Blair/DerCor_WGR_Final/trimmed_files

# Run FastQC on the raw reads
#fastqc $FQDIR/dc_${SAMPLE}_R1.fq.gz --outdir ./FastQC/raw_reads
#fastqc $FQDIR/dc_${SAMPLE}_R2.fq.gz --outdir ./FastQC/raw_reads

# Run Trimmomatic to trim for quality
#java -jar /share/pkg/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 -threads 16 $FQDIR/dc_${SAMPLE}_R1.fq.gz $FQDIR/dc_${SAMPLE}_R2.fq.gz \
#$TRIMDIR/paired/dc_${SAMPLE}_R1.paired.fq.gz $TRIMDIR/unpaired/dc_${SAMPLE}_R1.unpaired.fq.gz \
#$TRIMDIR/paired/dc_${SAMPLE}_R2.paired.fq.gz $TRIMDIR/unpaired/dc_${SAMPLE}_R2.unpaired.fq.gz \
#ILLUMINACLIP:/share/pkg/trimmomatic/0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:75

# Run FastQC on the trimmed reads
#fastqc $TRIMDIR/paired/dc_${SAMPLE}_R1.paired.fq.gz --outdir ./FastQC/trimmed_reads
#fastqc $TRIMDIR/paired/dc_${SAMPLE}_R2.paired.fq.gz --outdir ./FastQC/trimmed_reads

# Align trimmed reads to reference genome
#bwa mem -t 16 $REFDIR/$REF $TRIMDIR/paired/dc_${SAMPLE}_R1.paired.fq.gz $TRIMDIR/paired/dc_${SAMPLE}_R2.paired.fq.gz > $OUTDIR/dc_${SAMPLE}.sam
#rm $TRIMDIR/*/dc_${SAMPLE}*

# Conver to BAM and delete SAM file
#samtools view -S -b $OUTDIR/dc_${SAMPLE}.sam > $OUTDIR/dc_${SAMPLE}.bam
#rm $OUTDIR/dc_${SAMPLE}.sam

# Use SAMtools flagstat to get alignment statisctics:
#samtools flagstat $OUTDIR/dc_${SAMPLE}.bam > $METDIR/alignment_metrics/dc_${SAMPLE}_alignment_stats.txt

# Sort and index file, remove original BAM
#samtools sort -o $OUTDIR/dc_${SAMPLE}_SORT.bam $OUTDIR/dc_${SAMPLE}.bam
#samtools index $OUTDIR/dc_${SAMPLE}_SORT.bam
#rm $OUTDIR/dc_${SAMPLE}.bam

# Check mean depth
#source activate mosdepth
#cd ./metrics/depth_metrics/pre_dupe/
#mosdepth -n --fast-mode -Q 20 dc_${SAMPLE}_pre_dupe $OUTDIR/dc_${SAMPLE}_SORT.bam
#cd ../../../
#conda deactivate

# Remove duplicates
#java -Xms72G -Xmx90G -jar /share/pkg/picard/2.23.3/picard.jar MarkDuplicates \
#      I=$OUTDIR/dc_${SAMPLE}_SORT.bam \
#      O=$OUTDIR/dc_${SAMPLE}_SORT_DR.bam \
#      M=$METDIR/dupe_metrics/dc_${SAMPLE}_duplicate_metrics.txt \
#      ASSUME_SORT_ORDER=coordinate \
#      VALIDATION_STRINGENCY=SILENT \
#      REMOVE_DUPLICATES=TRUE \
#      MAX_RECORDS_IN_RAM=500000
      
#samtools index $OUTDIR/dc_${SAMPLE}_SORT_DR.bam

# Check mean depth and remove intermediate BAM
#source activate mosdepth
#cd ./metrics/depth_metrics/post_dupe/
#mosdepth -n --fast-mode -Q 20 dc_${SAMPLE}_post_dupe $OUTDIR/dc_${SAMPLE}_SORT_DR.bam
#cd ../../../
#conda deactivate

#rm $OUTDIR/dc_${SAMPLE}_SORT.bam*

# Add RG header
#java -jar -Xms72G -Xmx90G /share/pkg/picard/2.23.3/picard.jar AddOrReplaceReadGroups \
# I=$OUTDIR/dc_${SAMPLE}_SORT_DR.bam O=$OUTDIR/dc_${SAMPLE}_SORT_DR_RG.bam \
# SORT_ORDER=coordinate RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=dc_${SAMPLE} CREATE_INDEX=True
 
#samtools index $OUTDIR/dc_${SAMPLE}_SORT_DR_RG.bam

#rm $OUTDIR/dc_${SAMPLE}_SORT_DR.bam*

# Re-align around indels for use with ANGSD
#java -jar -Xms72G -Xmx90G /share/pkg/GATK/3.5/GenomeAnalysisTK.jar \
#  -T RealignerTargetCreator \
#  -R $REFDIR/$REF \
#  -I $OUTDIR/dc_${SAMPLE}_SORT_DR_RG.bam \
#  -o $OUTDIR/dc_${SAMPLE}.intervals \
#  -drf BadMate


#java -jar -Xms72G -Xmx90G /share/pkg/GATK/3.5/GenomeAnalysisTK.jar \
#  -T IndelRealigner \
#  -R $REFDIR/$REF \
#  -I $OUTDIR/dc_${SAMPLE}_SORT_DR_RG.bam \
#  -targetIntervals $OUTDIR/dc_${SAMPLE}.intervals \
#  --consensusDeterminationModel USE_READS  \
#  -o $OUTDIR/dc_${SAMPLE}_SORT_DR_RG_IR.bam
  
#rm $OUTDIR/dc_${SAMPLE}_SORT_DR_RG.bam*
#rm $OUTDIR/dc_${SAMPLE}.intervals

# Sort and index final time:
samtools sort $OUTDIR/dc_${SAMPLE}_SORT_DR_RG_IR.bam -o $OUTDIR/dc_${SAMPLE}_DR_RG_IR_SORT.bam
samtools index $OUTDIR/dc_${SAMPLE}_DR_RG_IR_SORT.bam