#!/bin/bash
#SBATCH -J DerCor_WGR_concat
#SBATCH -o ./logs/DerCor_WGR_concat_%A_%a.log
#SBATCH -e ./logs/DerCor_WGR_concat_%A_%a.err
#SBATCH -p cpu-long
#SBATCH -t 24:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=32000  # Requested Memory
#SBATCH --array=1

SAMPLE=$(sed -n ${SLURM_ARRAY_TASK_ID}p ./links/all_IDs.txt)
echo ${SAMPLE}
ls ./links/*${SAMPLE}*_1.fq.gz | wc -l

# R1
cat ./links/*${SAMPLE}*_1.fq.gz > ./concat/dc_${SAMPLE}_1.fq.gz

# R2
cat ./links/*${SAMPLE}*_2.fq.gz > ./concat/dc_${SAMPLE}_2.fq.gz

# R1
#cat /project/uma_lisa_komoroske/raw_reads/DerCor_WGR_all/DcWGR_HP/*${SAMPLE}*_1* \
#    /project/uma_lisa_komoroske/raw_reads/DerCor_WGR_all/DerCor_WGR_Round1/*${SAMPLE}*_1* \
#    /project/uma_lisa_komoroske/raw_reads/DerCor_WGR_all/DerCor_WGR_Round2/raw_data/*/*${SAMPLE}*_1* \
#    /project/uma_lisa_komoroske/raw_reads/DerCor_WGR_all/DerCor_WGR_Round3/usftp21.novogene.com/raw_data/*/*${SAMPLE}*_1* \
#    > ./concat_files/dc_${SAMPLE}_R1.fq.gz
    
# R2
#cat /project/uma_lisa_komoroske/raw_reads/DerCor_WGR_all/DcWGR_HP/*${SAMPLE}*_2* \
#    /project/uma_lisa_komoroske/raw_reads/DerCor_WGR_all/DerCor_WGR_Round1/*${SAMPLE}*_2* \
#    /project/uma_lisa_komoroske/raw_reads/DerCor_WGR_all/DerCor_WGR_Round2/raw_data/*/*${SAMPLE}*_2* \
#    /project/uma_lisa_komoroske/raw_reads/DerCor_WGR_all/DerCor_WGR_Round3/usftp21.novogene.com/raw_data/*/*${SAMPLE}*_2* \
#    > ./concat_files/dc_${SAMPLE}_R2.fq.gz