
# After making all the window assignments, I have 5 input files to feed to WZA
# To know what commands needed to be run, I used the following construction

Sinteractive
module load gcc
module load python
awk '{print FILENAME; nextfile}' *_WZAinput.csv | sed "s/input.csv//" | sed "s/2.//" | awk '{print "python general_WZA_script.py \
--correlations 2."$0"input.csv --summary_stat pvalue --window ID --MAF MAF --output 3."$0"output.csv --sep \",\""}'

# After getting 5 output files, I need to merge them

awk 'NR == 1 || FNR > 1' 3.* > 3.5_WZAoutputMerged.csv
