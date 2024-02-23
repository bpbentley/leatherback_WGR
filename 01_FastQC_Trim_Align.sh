#!/bin/bash
#SBATCH -J Duffy_Trim_Align
#SBATCH -o ./logs/align/Duffy_QC_Trim_Align_%a.log
#SBATCH -e ./logs/align/Duffy_QC_Trim_Align_%a.err
#SBATCH -p cpu-long
#SBATCH -t 96:00:00
#SBATCH -c 8  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=100000  # Requested Memory
#SBATCH --array=2 # Number of samples

module load bwa/0.7.17
module load samtools/1.9
module load python/3.11.0
module load fastqc/0.11.9
module load trimmomatic/0.39
module load picard/2.26.11
module load htslib/1.12
module load miniconda/22.11.1-1
module load java/11.0.2

#############################################################################################
### NOTE HERE THAT YOU NEED TO INSTALL GATK3 & MOSDEPTH THROUGH CONDA BEFORE RUNNING THIS ###
################ RUN THE NEXT COMMAND ONCE IN YOUR TERMINAL TO INSTALL: #####################
#############################################################################################

# conda create -n gatk3 -c conda-forge -c bioconda gatk
# Then gatk-register to the GenomeAnalysisTK.jar file (in the ./opt/ folder)
# While the conda environment is still active (conda activate gatk), install mosdepth:
# conda install mosdepth
# Should now be set up to use it

REFDIR=/nese/meclab/Shared/reference_genomes/ST_reference_genomes/rDerCor1_202105 # You need to download your reference genome to here
REF=rDerCor1.pri.cur.20210524.fasta
FQDIR=/nese/meclab/Blair/WGR_DerCor/data/duffy_data/concat # This is the directory where your raw data is: recommend using ln -s to the downloaded files
SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}'p' ./duffy_samples.txt) # This is the sample ID, have a text doc that is a list of the sample names
OUTDIR=/nese/meclab/Blair/WGR_DerCor/data/alignments/duffy # This is the folder your alignments will output to
METDIR=/nese/meclab/Blair/WGR_DerCor/metrics # This is where you save the metrics to
TRIMDIR=/nese/meclab/Blair/WGR_DerCor/trimmed # This is the output directory for trimmomatic (i.e. the files you'll use for input to BWA)

# Run FastQC on the raw reads
#fastqc $FQDIR/${SAMPLE}_R1.fq.gz --outdir ./FastQC/raw_reads
#fastqc $FQDIR/${SAMPLE}_R2.fq.gz --outdir ./FastQC/raw_reads

# Run Trimmomatic to trim for quality
java -jar /modules/apps/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 -threads 16 $FQDIR/${SAMPLE}_R1.fq.gz $FQDIR/${SAMPLE}_R2.fq.gz \
$TRIMDIR/paired/${SAMPLE}_R1.paired.fq.gz $TRIMDIR/unpaired/${SAMPLE}_R1.unpaired.fq.gz \
$TRIMDIR/paired/${SAMPLE}_R2.paired.fq.gz $TRIMDIR/unpaired/${SAMPLE}_R2.unpaired.fq.gz \
ILLUMINACLIP:/share/pkg/trimmomatic/0.39/adapters/NexteraPE-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:75 

# Run FastQC on the trimmed reads
fastqc $TRIMDIR/paired/${SAMPLE}_R1.paired.fq.gz --outdir ./FastQC/trimmed_reads
fastqc $TRIMDIR/paired/${SAMPLE}_R2.paired.fq.gz --outdir ./FastQC/trimmed_reads

# Align trimmed reads to reference genome
# Make sure reference has been indexed with BWA
bwa mem -t 16 $REFDIR/$REF $TRIMDIR/paired/${SAMPLE}_R1.paired.fq.gz $TRIMDIR/paired/${SAMPLE}_R2.paired.fq.gz > $OUTDIR/${SAMPLE}.sam
#rm $TRIMDIR/*/${SAMPLE}* # Might not want to run this - deletes the Trimmomatic outputs.

# Convert to BAM and delete SAM file
samtools view -S -b $OUTDIR/${SAMPLE}.sam > $OUTDIR/${SAMPLE}.bam
#rm $OUTDIR/${SAMPLE}.sam

# Use SAMtools flagstat to get alignment statisctics:
samtools flagstat $OUTDIR/${SAMPLE}.bam > $METDIR/alignment_metrics/${SAMPLE}_alignment_stats.txt

# Sort and index file, remove original BAM
samtools sort -o $OUTDIR/${SAMPLE}_SORT.bam $OUTDIR/${SAMPLE}.bam
samtools index $OUTDIR/${SAMPLE}_SORT.bam
#rm $OUTDIR/${SAMPLE}.bam

# Check mean depth
conda activate gatk3
cd ./metrics/depth_metrics/pre_dupe/ # Make sure this directory exists
mosdepth -n --fast-mode -Q 20 ${SAMPLE}_pre_dupe $OUTDIR/${SAMPLE}_SORT.bam
cd ../../../
conda deactivate

# Remove duplicates
java -Xms72G -Xmx90G -jar /modules/apps/picard/2.26.11/picard.jar MarkDuplicates \
      I=$OUTDIR/${SAMPLE}_SORT.bam \
      O=$OUTDIR/${SAMPLE}_SORT_DR.bam \
      M=$METDIR/dupe_metrics/${SAMPLE}_duplicate_metrics.txt \
      ASSUME_SORT_ORDER=coordinate \
      VALIDATION_STRINGENCY=SILENT \
      REMOVE_DUPLICATES=TRUE \
      MAX_RECORDS_IN_RAM=500000
      
samtools index $OUTDIR/${SAMPLE}_SORT_DR.bam

# Check mean depth and remove intermediate BAM
conda activate gatk3
cd ./metrics/depth_metrics/post_dupe/
mosdepth -n --fast-mode -Q 20 ${SAMPLE}_post_dupe $OUTDIR/${SAMPLE}_SORT_DR.bam
cd ../../../
conda deactivate

#rm $OUTDIR/${SAMPLE}_SORT.bam*

# Add RG header
java -jar -Xms72G -Xmx90G /modules/apps/picard/2.26.11/picard.jar AddOrReplaceReadGroups \
 I=$OUTDIR/${SAMPLE}_SORT_DR.bam O=$OUTDIR/${SAMPLE}_SORT_DR_RG.bam \
 SORT_ORDER=coordinate RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=${SAMPLE} CREATE_INDEX=True
 
samtools index $OUTDIR/${SAMPLE}_SORT_DR_RG.bam

#rm $OUTDIR/${SAMPLE}_SORT_DR.bam*

# Re-align around indels for use with ANGSD
# Need to check that this works and if there's a way to specify memory - I've usually used it with Java, but not sure how conda handles the memory requests
conda activate gatk3
gatk \
  -T RealignerTargetCreator \
  -R $REFDIR/$REF \
  -I $OUTDIR/${SAMPLE}_SORT_DR_RG.bam \
  -o $OUTDIR/${SAMPLE}.intervals \
  -drf BadMate


gatk \
  -T IndelRealigner \
  -R $REFDIR/$REF \
  -I $OUTDIR/${SAMPLE}_SORT_DR_RG.bam \
  -targetIntervals $OUTDIR/${SAMPLE}.intervals \
  --consensusDeterminationModel USE_READS  \
  -o $OUTDIR/${SAMPLE}_SORT_DR_RG_IR.bam

conda deactivate
  
#rm $OUTDIR/${SAMPLE}_SORT_DR_RG.bam*
#rm $OUTDIR/${SAMPLE}.intervals

# Sort and index final time:
samtools sort $OUTDIR/${SAMPLE}_SORT_DR_RG_IR.bam -o $OUTDIR/${SAMPLE}_DR_RG_IR_SORT.bam
samtools index $OUTDIR/${SAMPLE}_DR_RG_IR_SORT.bam