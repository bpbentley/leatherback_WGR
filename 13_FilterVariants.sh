#!/bin/bash
#SBATCH -J WGR_FilterVariants
#SBATCH -o ./logs/WGR_FilterVariants.log
#SBATCH -e ./logs/WGR_FilterVariants.err
#SBATCH -p cpu
#SBATCH -t 04:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=96000  # Requested Memory


module load vcftools/0.1.14

for CHR in {1..28}; do

input_vcf=/nese/meclab/Blair/WGR_DerCor/snpEff/tvcf/DerCor_SV_SUPER_${CHR}.vcf
fileprefix=$(echo SUPER_$CHR)

######################################
### Filter VCF for use with snpEff ###
######################################

### 
vcftools --vcf $input_vcf --min-meanDP 5 --max-meanDP 200 --recode --mac 1 --out /nese/meclab/Blair/WGR_DerCor/snpEff/fvcf/${fileprefix}

done