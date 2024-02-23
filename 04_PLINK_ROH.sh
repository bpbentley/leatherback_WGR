#!/bin/bash
#SBATCH -J ROH_PLINK_filt_DerCor
#SBATCH -o ./logs/ROH_PLINK_filt_DerCor.log
#SBATCH -e ./logs/ROH_PLINK_filt_DerCor.err
#SBATCH -p cpu
#SBATCH -t 04:00:00
#SBATCH -c 16  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=96000  # Requested Memory

####################
### Load modules ###
####################

module load plink/1.07

######################
### Set file paths ###
######################
INDIR=/nese/meclab/Blair/WGR_DerCor/ANGSD/for_PLINK/filt
PLINKDIR=/nese/meclab/Blair/WGR_DerCor/ROH/PLINK/filt


#####################################################
### Run PLINK ROH estimation for each chromosome: ###
#####################################################

for p in {1..28}; do
plink --tfile $INDIR/Filt_DerCor_SUPER_${p} --homozyg-snp 20 --homozyg-kb 50 \
 --homozyg-window-snp 20 --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.01 \
 --out $PLINKDIR/Filt_DerCor_SUPER_${p} --noweb
done

##########################
### Testing Paramaters ###
##########################

#for q in {0..100..10}; do
#for p in {1..28}; do
#plink --tfile $INDIR/'subset_DerCor_SUPER_'$p --homozyg-snp 20 --homozyg-kb 50 \
# --homozyg-window-snp 20 --homozyg-window-het 3 --homozyg-window-missing ${q} --homozyg-window-threshold 0.01 \
# --allow-extra-chr --out $PLINKDIR/test_MISS/MISS_${q}/subset_DerCor_${p}_MISS_${q}
#done
#done

# --homozyg-density 20 --homozyg-gap 1000

#for p in {1..28}; do
#for q in {10..200..10}; do
#plink --tfile $INDIR/'WGR_DerCorTop5_SUPER_'${p} --homozyg-snp ${q} --homozyg-kb 1 \
# --homozyg-window-snp ${q} --homozyg-window-het 1 --homozyg-window-missing 5 --homozyg-window-threshold 0.01 \
# --allow-extra-chr --out $PLINKDIR/test_SNP/SNP_${q}/DerCorTop5_SUPER_${p}_SNP_${q}
#done
#done