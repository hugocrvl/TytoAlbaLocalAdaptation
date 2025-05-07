#!/bin/bash

#SBATCH --time 04:00:00
#SBATCH --mem 20G
#SBATCH --cpus-per-task 10

module load gcc
module load bcftools
module load vcftools

bcftools concat -f output_Sphased.txt --output output_shapeit/Phased_shapeit.vcf

vcftools --vcf output_whatshap/Phased_whatshap.vcf --chr Super-Scaffold_2 --recode --out evaluation/Phased_whatshap_SC2.vcf
mv evaluation/Phased_whatshap_SC2.vcf.recode.vcf evaluation/Phased_whatshap_SC2.vcf
bcftools view -Oz --output evaluation/Phased_whatshap_SC2.vcf.gz evaluation/Phased_whatshap_SC2.vcf
bcftools index evaluation/Phased_whatshap_SC2.vcf.gz