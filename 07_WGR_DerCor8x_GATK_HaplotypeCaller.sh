#!/bin/bash
#SBATCH -J GATK_HapCaller
#SBATCH -o ./logs/GATK/GATK_HapCaller_%a.log
#SBATCH -e ./logs/GATK/GATK_HapCaller_%a.err
#SBATCH -p cpu-long
#SBATCH -t 48:00:00
#SBATCH -c 12  # Number of Cores per Task
#SBATCH --threads-per-core=2
#SBATCH --mem=64000  # Requested Memory
#SBATCH --array=1-40

# Load required modules
module load gatk/4.2.3.0
module load htslib/1.9
module load samtools/1.9
module load python/3.9.1

###############################################################################
# These values will need to be adjusted to fit your environment as well
# as on a per analysis run basis

# NOTE: Important!!!! Need to set the min and max depth for the genome ID in the python script filterVCF_20200723.py
	#   Good idea to try to modify this script to update the python script prior to running. 

# Adjust NAME to whatever prefix you want your output files to have
SAMPLE=$(printf $(sed -n ${SLURM_ARRAY_TASK_ID}'p' ./high_coverage_samples.txt)) # unique name needed if running multiple pipelines on different data, as files will be written to same directory, and will over-write if not unique.
NAME=${SAMPLE}

# Directories
BAMDIR=/nese/meclab/Blair/WGR_DerCor/data/alignments
OUTDIR=/nese/meclab/Blair/WGR_DerCor/GATK/GVCF
REFDIR=/nese/meclab/Shared/reference_genomes/ST_reference_genomes/rDerCor1_202105

# Files
REFERENCE=rDerCor1.pri.cur.20210524.fasta  # reference must be indexed with samtools faidx, and data dictionary (samtools dict) first
BAM=${NAME}_DR_RG_IR_SORT.bam # samtools index bam file first!

# Output Files
MYLOG=${OUTDIR}/${NAME}_gatk.log
GVCF=${OUTDIR}/${NAME}.g.vcf.gz

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
    -I ${BAMDIR}/${BAM} \
    -O ${GVCF} \
    >> ${MYLOG} 2>&1
