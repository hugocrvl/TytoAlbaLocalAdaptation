#!/bin/bash
#SBATCH --cpus-per-task 10
#SBATCH --mem 20G
#SBATCH --time 02:00:00

module load gcc
module load vcftools
module load bcftools
module load r

vcftools --vcf Complete_EU_Geo_filter_3.vcf --exclude-bed all_mask_150_90.bed --keep ID_kept.txt --max-missing 0.95 \
 	--not-chr Super-Scaffold_200000105744 --not-chr Super-Scaffold_2000008 \
        --not-chr Super-Scaffold_20000014 --not-chr Super-Scaffold_100000100381 --not-chr Super-Scaffold_200000100656 \
        --not-chr Super-Scaffold_15 --not-chr Super-Scaffold_200000178 --not-chr Super-Scaffold_24 \
        --not-chr Super-Scaffold_100000100064 --not-chr Super-Scaffold_37 \
        --not-chr Super-Scaffold_100000100052 --not-chr Super-Scaffold_100000100078 --not-chr Super-Scaffold_13 \
        --not-chr Super-Scaffold_4 --not-chr Super-Scaffold_25 --not-chr Super-Scaffold_42 \
	--min-alleles 2 --max-alleles 2 --recode --out phasing/Complete_EU_filters_masked_woMAF.vcf
mv phasing/Complete_EU_filters_masked_woMAF.vcf.recode.vcf phasing/Complete_EU_filters_masked_woMAF.vcf

bcftools convert -O b phasing/Complete_EU_filters_masked_woMAF.vcf > phasing/Complete_EU_filters_masked_woMAF.bcf
bcftools index phasing/Complete_EU_filters_masked_woMAF.bcf

