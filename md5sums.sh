#!/bin/bash
#BSUB -J "DerCor_WGR_align[1-167]"
#BSUB -o ./logs/md5sums/md5sum_%I.log
#BSUB -e ./logs/md5sums/md5sum_%I.err
#BSUB -q long
#BSUB -W 24:00
#BSUB -R rusage[mem=4000]
#BSUB -R span[hosts=1]
#BSUB -n 16

SAMPLE=$(sed -n ${LSB_JOBINDEX}p ./all_bam_files.txt)

md5sum ${SAMPLE} >> md5sums_ghpcc.txt