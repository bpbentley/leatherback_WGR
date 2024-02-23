#!/bin/bash
#SBATCH -J NGSrelate
#SBATCH -o ./logs/NGSrelate/Run_NGSrelate.log
#SBATCH -e ./logs/NGSrelate/Run_NGSrelate.err
#SBATCH -p cpu-long
#SBATCH -t 48:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=96000  # Requested Memory

module load angsd/0.935

INDIR=/nese/meclab/Blair/WGR_DerCor/ANGSD/NGSrelate
OUTDIR=/nese/meclab/Blair/WGR_DerCor/ngsRelate
CHR=SUPER_${SLURM_ARRAY_TASK_ID}
REFDIR=/nese/meclab/Shared/reference_genomes/ST_reference_genomes/rDerCor1_202105
REF=rDerCor1.pri.cur.20210524.fasta
NGSREL=/nese/meclab/Shared/bin/ngsRelate

#angsd -out $OUTDIR/DerCor_NGSrelate_${CHR} -b ./filtered_alignments.txt -gl 1 -domajorminor 1 -snp_pval 1e-6 -domaf 1 -minmaf 0.05 -doGlf 3 -nThreads 16 -r $CHR -ref $REFDIR/$REF

# Combine first
zcat $INDIR/DerCor_NGSrelate.mafs.gz | cut -f6 |sed 1d > $INDIR/freq

$NGSREL/ngsRelate  -g $INDIR/DerCor_NGSrelate.glf.gz -n 100 -f $INDIR/freq  -O $OUTDIR/newres