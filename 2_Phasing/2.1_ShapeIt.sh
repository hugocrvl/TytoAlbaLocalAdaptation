#!/bin/bash

#SBATCH --time 03:00:00
#SBATCH --cpus-per-task 10
#SBATCH --mem 20G
#SBATCH --array=0-54

module load gcc
module load vcftools
module load bcftools
module load shapeit4

while read scaffold; do scaffold_list+=($scaffold); done < '../scaffolds.txt'
while read output_shapeit; do output_shapeit_list+=($output_shapeit); done < 'output_Sphased.txt'

shapeit4.2 --input output_whatshap/Phased_whatshap.vcf.gz --region ${scaffold_list[$SLURM_ARRAY_TASK_ID]} --output ${output_shapeit_list[$SLURM_ARRAY_TASK_ID]}
