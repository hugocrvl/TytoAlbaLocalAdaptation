Scripts necessary to phase genetic variants.

1. `1.1_WhatsHap.sh` phases per individual and `1.2_MergingIndividualPhasing.sh` merges the output. Additionnally, the latter script filters for --maf 0.05 with `VCFTools`. 
2. `2.1_ShapeIt.sh` phases per scaffold and `2.2_MergingScaffoldPhaing.sh` merges the output.
3. `3_Evaluation` directory contains one script to assess the Switch Error Rate along the Super-Scaffold 2 (longest of the current assembly).

The files in the `0_AdditionalFiles` directory were used during the phasing process as input and output variables

**Disclaimer**: The directory names referenced in the scripts may not exactly match the directory names used in this GitHub repository.
