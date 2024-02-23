#!/bin/bash
#BSUB -J "WGR_DerCorTop5_GATK_JointGeno[1-28]"
#BSUB -e ./logs/07_GATK_PLINK/rDerCorTop5_GATK_JointGeno.%I.err
#BSUB -o ./logs/07_GATK_PLINK/rDerCorTop5_GATK_JointGeno.%I.log
#BSUB -W 72:00
#BSUB -n 16
#BSUB -R rusage[mem=6000]
#BSUB -R span[hosts=1]
#BSUB -q long

# Load required modules
module load gatk/4.1.8.1
module load htslib/1.9
module load samtools/1.4.1
module load python3/3.5.0
module load htslib/1.9

INDIR=/project/uma_lisa_komoroske/Blair/WGR_DerCor/GATK/GVCF
DERCOR=/project/uma_lisa_komoroske/Blair/rDerCor1/GATK_out/genome/temp
OUTDIR=/project/uma_lisa_komoroske/Blair/WGR_DerCor/GATK/VCF
REFDIR=/project/uma_lisa_komoroske/Blair/refs/rDerCor1_20210524
REF=rDerCor1.pri.cur.20210524.fasta
SCAFFOLDLIST=rDerCor1.pri.cur.20210524.scaffold.list.txt

NUM=$(printf %02d ${LSB_JOBINDEX})
CHR=$(head -n ${NUM} ${REFDIR}/${SCAFFOLDLIST} | tail -n 1)
CHR=$(echo ${CHR} | awk -F " " '{ print $1 }')
echo ${CHR}

source /share/pkg/condas/2018-05-11/bin/activate && conda activate gatk_4.1.8.1

### Compile the gVCFs into a single database (this replaces the CombineGVCF function)
#### Can still use CombineGVCF but this is quicker and more efficient as per the best practices
gatk \
     GenomicsDBImport \
    --java-options "-Xms60G -Xmx90G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -V $INDIR/s_100028.g.vcf.gz \
    -V $INDIR/s_11156.g.vcf.gz \
    -V $INDIR/s_14559.g.vcf.gz \
    -V $INDIR/s_14560.g.vcf.gz \
    -V $INDIR/s_37134.g.vcf.gz \
    -V $DERCOR/rDerCor_genome.g.vcf.gz \
    --genomicsdb-workspace-path ${INDIR}/DerCorTop5DB_${CHR} \
    -L ${CHR}
   
### Perform joint-genotyping of the individuals
#### Equivalent to running all the samples through the ANGSD pipeline and generating a SNP-list
gatk \
    GenotypeGVCFs \
    --java-options "-Xms60G -Xmx90G -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" \
    -R $REFDIR/$REF \
    -V gendb://${INDIR}/DerCorTop5DB_${CHR} \
    -O ${OUTDIR}/DerCor_Top5_${CHR}.vcf \
    -L ${CHR}
    
conda deactivate && conda deactivate
