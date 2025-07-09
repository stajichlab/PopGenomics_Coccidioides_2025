Coccidioides Population Genomics Analysis 
===

Subfolders organized around mapping to a specific reference and species sample set.

1. Dataset
   * [_C. immitis_ strains](Genotyping/C_immitis.samples.csv)
   * [_C. posadasii_ strains](C_posadasii.samples.csv)
   * [all Coccidioides strains](Coccidioides.samples.csv) 
2. Primary Genotyping and Comparison Analyses
   * Variant calling (reads to VCF)
     * [_C. immitis_ RS](C_immitis_ref_RS) _C. immitis_ only aligned to RS genome
     * [_C. posadasii_ str Silveira 2022](C_posadasii_refSilv2022) _C. posadasii_ only aligned to Silveira ([2022 Nanopore genome](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_018416015.2/))
     * [All strains to C. immitis RS](all_C_immitis_RS] - all Coccidioides aligned to RS
     * [All strains to C. immitis WA 211](all_C_immitis_ref_WA211] - all Coccidioides aligned to Washington 211 asm (this Ci. strain does not have evidence of introgressions)
     * [All strains to C. posadasii str Silveira 2022](all_C_posadasii_refSilv2022] - all Coccidioides aligned to Silveira 2022 asm
   * CNV
     * [Ci RS ref CNV_plots](all_C_immitis_ref_RS/results/CNV_plots/]
     * [Ci WA211 ref CNV_plots](all_C_immitis_ref_WA211/results/CNV_plots/]
     * [Cp Silv2022 ref CNV_plots](all_C_posadasii_refSilv2022/results/CNV_plots]
3. Population and Genomewide analyses
  * Population Structure, Introgression testing
  * Scans for selection and summary statistics
  * Identifying outlier genes with high substitution counts
  * GWAS tests for association with strain traits
  * Drug resistance alleles
