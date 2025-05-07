#!/bin/bash

#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --job-name=EvalPhase
#SBATCH --array=0-73

module load gcc
module load vcftools
module load bcftools
module load whatshap
module load shapeit4/4.1.3

while read name; do names+=($name); done < '../ID.txt'

vcftools --vcf Phased_whatshap_SC2.vcf \
		 --out IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}only \
		 --recode \
		 --indv ${names[$SLURM_ARRAY_TASK_ID]}

vcftools --vcf Phased_whatshap_SC2.vcf \
		 --out IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}free \
		 --recode \
		 --remove-indv ${names[$SLURM_ARRAY_TASK_ID]}


mv IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}only.recode.vcf IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}only.vcf
mv IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}free.recode.vcf IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}free.vcf
# rm *log

whatshap unphase IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}only.vcf > IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}only.unphased.vcf

bcftools view -Oz --output IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}free.vcf.gz IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}free.vcf
bcftools index IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}free.vcf.gz
bcftools view -Oz --output IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}only.unphased.vcf.gz IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}only.unphased.vcf
bcftools index IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}only.unphased.vcf.gz

bcftools query -s ${names[$SLURM_ARRAY_TASK_ID]} -f '[%PS\t]\n' Phased_whatshap_SC2.vcf | awk '{print $1}' | sed 's/\./NA/g' > IndivPS/${names[$SLURM_ARRAY_TASK_ID]}.PS

bcftools merge IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}free.vcf.gz IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}only.unphased.vcf.gz --output IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}unphased.vcf

bcftools view -Oz --output IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}unphased.vcf.gz IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}unphased.vcf
bcftools index IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}unphased.vcf.gz

# Phase with shapeit

shapeit4 --input IndVCF/AllInds_WhatsapPhased_mis95_Super_Scaffold_2_${names[$SLURM_ARRAY_TASK_ID]}unphased.vcf.gz \
-R Super-Scaffold_2 \
--output Out_ShapeIt/AllInds_ShapeItPhased_${names[$SLURM_ARRAY_TASK_ID]}.vcf

vcftools --vcf Out_ShapeIt/AllInds_ShapeItPhased_${names[$SLURM_ARRAY_TASK_ID]}.vcf \
		 --out Out_ShapeIt/${names[$SLURM_ARRAY_TASK_ID]}_ShapeItPhased \
		 --recode \
		 --indv ${names[$SLURM_ARRAY_TASK_ID]}

mv Out_ShapeIt/${names[$SLURM_ARRAY_TASK_ID]}_ShapeItPhased.recode.vcf Out_ShapeIt/${names[$SLURM_ARRAY_TASK_ID]}_ShapeItPhased.vcf

bcftools view -Oz --output Out_ShapeIt/${names[$SLURM_ARRAY_TASK_ID]}_ShapeItPhased.vcf.gz Out_ShapeIt/${names[$SLURM_ARRAY_TASK_ID]}_ShapeItPhased.vcf
bcftools index Out_ShapeIt/${names[$SLURM_ARRAY_TASK_ID]}_ShapeItPhased.vcf.gz

# Run SwithError

switchError/bin/swithError --gen Phased_whatshap_SC2.vcf.gz --hap Out_ShapeIt/${names[$SLURM_ARRAY_TASK_ID]}_ShapeItPhased.vcf.gz --ps IndivPS/${names[$SLURM_ARRAY_TASK_ID]}.PS --out OutSwith/${names[$SLURM_ARRAY_TASK_ID]} --reg Super-Scaffold_2

rm IndVCF/*${names[$SLURM_ARRAY_TASK_ID]}*
