#!/bin/bash
#SBATCH -J concat_Duffy_data
#SBATCH -o ./logs/concat_Duffy_data.log
#SBATCH -e ./logs/concat_Duffy_data.err
#SBATCH -p gpu-long
#SBATCH -t 12:00:00
#SBATCH -c 4  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=16000  # Requested Memory

cat ./data/duffy_data/DcBR3ECB*R1* > ./data/duffy_data/concat/DcBR3ECB_R1.fq.gz
cat ./data/duffy_data/DcBR3ECB*R2* > ./data/duffy_data/concat/DcBR3ECB_R2.fq.gz

cat ./data/duffy_data/LeatherHatchDNA*R1* > ./data/duffy_data/concat/LeatherHatchDNA_R1.fq.gz
cat ./data/duffy_data/LeatherHatchDNA*R2* > ./data/duffy_data/concat/LeatherHatchDNA_R2.fq.gz