#!/bin/bash

#SBATCH --time 24:00:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 20
#SBATCH --mem 50G
#SBATCH --array=0-73

module load gcc
module load whatshap
module load vcftools
module load bcftools

while read id; do id_list+=($id); done < 'ID.txt'
while read bam; do bam_list+=($bam); done < 'bam_files.txt'
while read output_all; do output_all_list+=($output_all); done < 'output_Wphased_all.txt'
while read output_ind; do output_ind_list+=($output_ind); done < 'output_Wphased_ind.txt'

whatshap phase -o ${output_all_list[$SLURM_ARRAY_TASK_ID]} --no-reference --sample ${id_list[$SLURM_ARRAY_TASK_ID]} Complete_EU_filters_masked_woMAF.vcf ${bam_list[$SLURM_ARRAY_TASK_ID]}

vcftools --vcf ${output_all_list[$SLURM_ARRAY_TASK_ID]} --indv ${id_list[$SLURM_ARRAY_TASK_ID]} --recode --out ${output_ind_list[$SLURM_ARRAY_TASK_ID]}
mv ${output_ind_list[$SLURM_ARRAY_TASK_ID]}.recode.vcf ${output_ind_list[$SLURM_ARRAY_TASK_ID]}
rm ${output_all_list[$SLURM_ARRAY_TASK_ID]}

bcftools view -Oz --output ${output_ind_list[$SLURM_ARRAY_TASK_ID]}.gz ${output_ind_list[$SLURM_ARRAY_TASK_ID]}
bcftools index ${output_ind_list[$SLURM_ARRAY_TASK_ID]}.gz
