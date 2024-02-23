#!/bin/bash
#SBATCH -J DerCor_WGR_Paired_SFS
#SBATCH -o ./logs/SAF_gen/DerCor_Paired_SFS_%a_%A.log
#SBATCH -e ./logs/SAF_gen/DerCor_Paired_SFS_%a_%A.err
#SBATCH -p cpu-long
#SBATCH -t 48:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=164000  # Requested Memory
#SBATCH --array=11

POP=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./popmaps/Populations.txt)
OUTDIR=/nese/meclab/Blair/WGR_DerCor/ANGSD/SAF
OUTDIR2=/nese/meclab/Blair/WGR_DerCor/ANGSD/SFS/2DSFS
REFDIR=/nese/meclab/Shared/reference_genomes/ST_reference_genomes/rDerCor1_202105
REF=rDerCor1.pri.cur.20210524.fasta

module load angsd/0.935

echo $POP

#angsd -b ./popmaps/$POP'_bams.txt' -ref ${REFDIR}/$REF -anc ${REFDIR}/$REF -out ${OUTDIR}/$POP \
#      -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
#      -minMapQ 20 -minQ 20 -minInd 2 -setMinDepth 7 -setMaxDepth 30 -doCounts 1 \
#      -GL 1 -doSaf 1
      
#realSFS ${OUTDIR}/$POP.saf.idx -fold 1 > ${OUTDIR2}/$POP.folded.sfs

for NUM in {1..12}; do
POP2=$(sed -n ${NUM}p ./popmaps/Populations.txt)

realSFS ${OUTDIR}/$POP.saf.idx ${OUTDIR}/$POP2.saf.idx -fold 1 > ${OUTDIR2}/$POP.$POP2.folded.sfs

done