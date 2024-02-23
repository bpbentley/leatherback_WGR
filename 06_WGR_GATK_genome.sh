#!/bin/bash
#SBATCH -J GATK_genome
#SBATCH -o ./logs/GATK/GATK_genome_%a_%j.log
#SBATCH -e ./logs/GATK/GATK_genome_%a_%j.err
#SBATCH -p cpu-long
#SBATCH -t 72:00:00
#SBATCH -c 12  # Number of Cores per Task
#SBATCH --threads-per-core=2
#SBATCH --mem=96000  # Requested Memory
#SBATCH --array=1-28

# Load required modules
module load gatk/4.2.3.0
module load htslib/1.9
module load samtools/1.9
module load python/3.9.1
module load miniconda/4.11.0

###############################################################################
# These values will need to be adjusted to fit your environment as well
# as on a per analysis run basis

# NOTE: Important!!!! Need to set the min and max depth for the genome ID in the python script filterVCF_20200723.py
	#   Good idea to try to modify this script to update the python script prior to running. 

# Adjust NAME to whatever prefix you want your output files to have
NAME=$(sed -n 1p ./high_coverage_samples.txt | cut -f1) # unique name needed if running multiple pipelines on different data, as files will be written to same directory, and will over-write if not unique.

# Directories
SCRIPTDIR=/nese/meclab/Blair/WGR_DerCor/scripts/gatk_in
BAMDIR=/nese/meclab/Blair/WGR_DerCor/data/alignments
OUTDIR=/nese/meclab/Blair/WGR_DerCor/GATK/${NAME}
TEMPDIR=/nese/meclab/Blair/WGR_DerCor/GATK/${NAME}/temp
REFDIR=/nese/meclab/Shared/reference_genomes/ST_reference_genomes/rDerCor1_202105

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
NUM=$(printf %02d ${SLURM_ARRAY_TASK_ID})
CHR=$(head -n ${NUM} ${REFDIR}/${SCAFFOLDLIST} | tail -n 1)
# NC_045784.1 Phocoena sinus isolate mPhoSin1 chromosome X, mPhoSin1.pri, whole genome shotgun sequence
CHR=$(echo ${CHR} | awk -F " " '{ print $1 }')
echo ${CHR}

# Output Files
MYLOG=${OUTDIR}/${NAME}_gatk.${NUM}.log
GVCF=${TEMPDIR}/${NAME}_${NUM}.g.vcf.gz
VCF=${TEMPDIR}/${NAME}_${NUM}.vcf.gz
TVCF=${TEMPDIR}/${NAME}_${NUM}.TrimAlt.vcf.gz
AVCF=${TEMPDIR}/${NAME}_${NUM}.AddAnnot.vcf.gz
FVCF=${OUTDIR}/${NAME}_${NUM}.Filter.vcf.gz

###############################################################################

# Start pipeline
# HaplotypeCaller
echo -e "[$(date "+%Y-%m-%d %T")] Starting HaplotypeCaller" >> ${MYLOG}
gatk \
    HaplotypeCaller \
    --java-options "-Xms60G -Xmx90G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -R ${REFDIR}/${REFERENCE} \
    -ERC BP_RESOLUTION \
    -mbq 20 \
    -L ${CHR} \
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
    --java-options "-Xms60G -Xmx90G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -R ${REFDIR}/${REFERENCE} \
    --include-non-variant-sites \
    --standard-min-confidence-threshold-for-calling 0 \
    -L ${CHR} \
    -V ${GVCF} \
    -O ${VCF} \
    >> ${MYLOG} 2>&1

#############
# SelectVariants
# --remove-unused-alternates replaces -trimAlternates
echo -e "[$(date "+%Y-%m-%d %T")] Starting SelectVariants" >> ${MYLOG}
gatk \
    SelectVariants \
    --java-options "-Xms60G -Xmx90G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -R ${REFDIR}/${REFERENCE} \
    --remove-unused-alternates \
    -L ${CHR} \
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
#     --java-options "-Xms60G -Xmx90G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
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
echo -e "[$(date "+%Y-%m-%d %T")] Starting Python Script" >> ${MYLOG}
python3 ${PSCRIPT} ${TVCF} ${MINDEPTH} ${MAXDEPTH} | bgzip > ${FVCF}
tabix -p vcf ${FVCF}


echo -e "[$(date "+%Y-%m-%d %T")] Finished Pipeline" >> ${MYLOG}
wait

##############################################################################
# WinHet script
# uses python, and module "pysam" for reading, manipulating and writing genomic data sets.
# to install pysam, source your python venv
conda activate turtles
# pip install pysam

CHR=$(head -n ${NUM} ${CHRLENGTHS} | tail -n 1)

python ${WSCRIPT} ${FVCF} ${WINSIZE} ${STEPSIZE} ${CHR} ${CHRLENGTHS}

conda deactivate