#!/bin/bash

#SBATCH --time 02:00:00
#SBATCH --cpus-per-task 30
#SBATCH --mem 30G

module load gcc
module load bcftools
module load vcftools

echo Merging has started

bcftools merge -Oz -l output_Wphased_ind_compress.txt --output output_whatshap/Phased_whatshap.vcf

vcftools --vcf output_whatshap/Phased_whatshap.vcf --maf 0.05 --recode --out output_whatshap/Phased_whatshap.vcf
mv output_whatshap/Phased_whatshap.vcf.recode.vcf output_whatshap/Phased_whatshap.vcf

bcftools view -Oz --output output_whatshap/Phased_whatshap.vcf.gz output_whatshap/Phased_whatshap.vcf
bcftools index output_whatshap/Phased_whatshap.vcf.gz
