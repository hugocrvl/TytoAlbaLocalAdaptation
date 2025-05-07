# First I needed to compute the MAF for each SNP
# Therefore I computed the allelic frequencies for each marker using VCFTools
# I use the phased VCF as I computed the RDA on this one

module load gcc
module load vcftools
vcftools --vcf ../Complete_EU_masked_phased.vcf --freq --out 1.1_AlleleFreq

# Then for each SNP I compute the MAF:
# 1. I remove the first line
# 2. For the 5th element of each line
# 3. I select the second element after the ':' character
# 4. If this element is higher than 0.5, I calculate 1-x, otherwise I print x directly
# 5. I redirect everything in the 1.2_MAF.txt file

awk '(NR>1)' 1.1_AlleleFreq.frq | awk '{print($5)}'| cut -d':' -f2 | awk '{if($1>0.5) print 1-$1; else print $1}' > 1.2_MAF.txt

# After making all the window assignments, I have 5 input files to feed to WZA
# To know what commands needed to be run, I used the following construction

Sinteractive
module load gcc
module load python
awk '{print FILENAME; nextfile}' *_WZAinput.csv | sed "s/input.csv//" | sed "s/2.//" | awk '{print "python general_WZA_script.py \
--correlations 2."$0"input.csv --summary_stat pvalue --window ID --MAF MAF --output 3."$0"output.csv --sep \",\""}'

# After getting 5 output files, I need to merge them but they all have the same header
# I used this command line but I don't have any idea why it worked 

awk 'NR == 1 || FNR > 1' 3.* > 3.5_WZAoutputMerged.csv
