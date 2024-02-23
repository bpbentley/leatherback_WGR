#!/bin/bash
#SBATCH -J wpac_concat_PCangsd_in
#SBATCH -o ./logs/wpac_concat_PCangsd_in_%a.log
#SBATCH -e ./logs/wpac_concat_PCangsd_in_%a.err
#SBATCH -p cpu-long
#SBATCH -t 12:00:00
#SBATCH -c 12  # Number of Cores per Task
#SBATCH --threads-per-core=1
#SBATCH --mem=32000  # Requested Memory

zcat /nese/meclab/Blair/WGR_DerCor/ANGSD/for_PCAngsd/wpac/DerCor_wpac_SUPER_1.beagle.gz > /nese/meclab/Blair/WGR_DerCor/ANGSD/for_PCAngsd/wpac/DerCor_wpac.beagle

for q in {2..28}; do
zcat /nese/meclab/Blair/WGR_DerCor/ANGSD/for_PCAngsd/DerCor_wpac_SUPER_${q}.beagle.gz | sed '1d' >> /nese/meclab/Blair/WGR_DerCor/ANGSD/for_PCAngsd/wpac/DerCor_wpac.beagle
done


gzip /nese/meclab/Blair/WGR_DerCor/ANGSD/for_PCAngsd/wpac/DerCor_wpac.beagle