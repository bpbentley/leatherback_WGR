#!/bin/bash
#SBATCH -J rDerCor_GATK_exons
#SBATCH -o ./logs/GATK_exons/GATK_exons_%a.log
#SBATCH -e ./logs/GATK_exons/GATK_exons_%a.err
#SBATCH -p cpu-long
#SBATCH -t 72:00:00
#SBATCH -c 12  # Number of Cores per Task
#SBATCH --mem=20000  # Requested Memory
#SBATCH --array=2-40 #(N=40)

# Load required modules
module load gatk/4.2.3.0
module load htslib/1.12
module load samtools/1.13+py3.8.12

###############################################################################
# These values will need to be adjusted to fit your environment as well
# as on a per analysis run basis

# NOTE: Important!!!! Need to set the min and max depth for the genome ID in the python script filterVCF_20200723.py
	#   Good idea to try to modify this script to update the python script prior to running. 

# Adjust NAME to whatever prefix you want your output files to have
NAME=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./high_coverage_samples.txt | cut -f1) # unique name needed if running multiple pipelines on different data, as files will be written to same directory, and will over-write if not unique.

# Directories
SCRIPTDIR=/nese/meclab/Blair/WGR_DerCor/scripts/gatk_in
BAMDIR=/nese/meclab/Blair/WGR_DerCor/data/alignments
#mkdir /project/uma_lisa_komoroske/Blair/WGR_for_genome/GATK_out/CheMyd/${NAME}/exons
#mkdir /project/uma_lisa_komoroske/Blair/WGR_for_genome/GATK_out/CheMyd/${NAME}/exons/temp
OUTDIR=/nese/meclab/Blair/WGR_DerCor/GATK/exons
TEMPDIR=/nese/meclab/Blair/WGR_DerCor/GATK/exons
REFDIR=/nese/meclab/Shared/reference_genomes/ST_reference_genomes/rDerCor1_202105
#WINDIR=/project/uma_lisa_komoroske/Blair/refs/rCheMyd1_20210528/exons

# Files
REFERENCE=rDerCor1.pri.cur.20210524.fasta  # reference must be indexed with samtools faidx, and data dictionary (samtools dict) first
SCAFFOLDLIST=rDerCor1.pri.cur.20210524.fasta.scaffold.list.txt # list must have ">" removed from before each scaffold! SEE README FILE FOR HOW TO GENERATE THE SCAFFOLDLIST FROM THE FASTA FILE
CHRLENGTHS=${REFDIR}/rDerCor1.pri.cur.20210524.fasta.scaffold.lengths.txt # SEE README FILE FOR HOW TO GENERATE THE CHROMOSOME LENGTHS FILE FROM THE FASTA INDEX FILE
BAM=${NAME}_DR_RG_IR_SORT.bam # samtools index bam file first!

# scripts
PSCRIPT=${SCRIPTDIR}/filterVCF_010920.py # python script for filtering vcf file
WSCRIPT=${SCRIPTDIR}/slidingWindowHet_010920.py # python script for creating sliding windows and counting heterozygotes per window. Specify window size and step size:

# values
WINSIZE=100000
STEPSIZE=100000
# minimum and maximum depth of coverage in BAM file (usually set to 1/3x and 2x average depth of coverage; assumes one sample per BAM file). 
MINDEPTH=$(sed -n 1p ./high_coverage_samples.txt | cut -f3)
MAXDEPTH=$(sed -n 1p ./high_coverage_samples.txt | cut -f4)

###############################################################################
# Get working scaffold based on array number
NUM=$(printf %02d ${LSB_JOBINDEX})
CHR=$(head -n ${NUM} ${REFDIR}/${SCAFFOLDLIST} | tail -n 1)
# NC_045784.1 Phocoena sinus isolate mPhoSin1 chromosome X, mPhoSin1.pri, whole genome shotgun sequence
CHR=$(echo ${CHR} | awk -F " " '{ print $1 }')
echo ${CHR}

# Output Files
MYLOG=${OUTDIR}/${NAME}_gatk.log
GVCF=${TEMPDIR}/${NAME}.g.vcf.gz
VCF=${TEMPDIR}/${NAME}.vcf.gz
TVCF=${TEMPDIR}/${NAME}.TrimAlt.vcf.gz
AVCF=${TEMPDIR}/${NAME}.AddAnnot.vcf.gz
FVCF=${OUTDIR}/${NAME}.Filter.vcf.gz

###############################################################################

# Start pipeline
# HaplotypeCaller
echo -e "[$(date "+%Y-%m-%d %T")] Starting HaplotypeCaller" >> ${MYLOG}
gatk \
    HaplotypeCaller \
    --java-options "-Xmx80G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -R ${REFDIR}/${REFERENCE} \
    -ERC BP_RESOLUTION \
    -mbq 20 \
    -L ${REFDIR}/rDerCor1.20210524.exons.renamed.trim.bed \
    -I ${BAMDIR}/${BAM} \
    -O ${GVCF} \
    --output-mode EMIT_ALL_ACTIVE_SITES \
    >> ${MYLOG} 2>&1

#############
# GenotypeGVCFs
# GATK3->GATK4 Update
# --include-non-variant-sites replaces -allSites
# --standard-min-confidence-threshold-for-calling replaces -stand_call_conf
echo -e "[$(date "+%Y-%m-%d %T")] Starting GenotypeGVCFs" >> ${MYLOG}
gatk \
    GenotypeGVCFs \
    --java-options "-Xmx80G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -R ${REFDIR}/${REFERENCE} \
    --include-non-variant-sites \
    --standard-min-confidence-threshold-for-calling 0 \
    -L ${REFDIR}/rDerCor1.20210524.exons.renamed.trim.bed \
    -V ${GVCF} \
    -O ${VCF} \
    >> ${MYLOG} 2>&1

#############
# SelectVariants
# --remove-unused-alternates replaces -trimAlternates
echo -e "[$(date "+%Y-%m-%d %T")] Starting SelectVariants" >> ${MYLOG}
gatk \
    SelectVariants \
    --java-options "-Xmx80G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -R ${REFDIR}/${REFERENCE} \
    --remove-unused-alternates \
    -L ${REFDIR}/rDerCor1.20210524.exons.renamed.trim.bed \
    -V ${VCF} \
    -O ${TVCF} \
    >> ${MYLOG} 2>&1

#############
# VariantAnnotator
# VariantType removed from GATK4 per
# https://gatkforums.broadinstitute.org/gatk/discussion/13500/a-varianttype-annotation-still-available-in-gatk4-haplotypecaller
# This is not functioning properly, and ends up changing the coverage per allele to 0,0 for all heterozygotes, so that there are none that pass filter in the next step.  Skip for now and run the final filter on the TrimAlt.vcf.gz files

# echo -e "[$(date "+%Y-%m-%d %T")] Starting VariantAnnotator" >> ${MYLOG}
# gatk \
#     VariantAnnotator \
#     --java-options "-Xmx80G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
#     -R ${REFDIR}/${REFERENCE} \
#     -G StandardAnnotation \
#     -L ${CHR} \
#     -V ${TVCF} \
#     -O ${AVCF} \
#     >> ${MYLOG} 2>&1

# python script to filter VCF (from AddAnnot.vcf.gz files)
# echo -e "[$(date "+%Y-%m-%d %T")] Starting Python Script" >> ${MYLOG}
# python3 ${PSCRIPT} ${AVCF} | bgzip > ${FVCF}
# tabix -p vcf ${FVCF}

#############
# Filter VCF
#python script to filter VCF (from TrimAlt.vcf.gz files)
#echo -e "[$(date "+%Y-%m-%d %T")] Starting Python Script" >> ${MYLOG}
#python3 ${PSCRIPT} ${TVCF} ${MINDEPTH} ${MAXDEPTH} | bgzip > ${FVCF}
#tabix -p vcf ${FVCF}


#echo -e "[$(date "+%Y-%m-%d %T")] Finished Pipeline" >> ${MYLOG}
#wait

##############################################################################
# WinHet script
# uses python, and module "pysam" for reading, manipulating and writing genomic data sets.
# to install pysam, source your python venv
#source /home/bb89a/python3/bin/activate
# pip install pysam

#CHR=$(head -n ${NUM} ${CHRLENGTHS} | tail -n 1)

#python ${WSCRIPT} ${FVCF} ${WINSIZE} ${STEPSIZE} ${CHR} ${CHRLENGTHS}

#deactivate
#conda deactivate && conda deactivate