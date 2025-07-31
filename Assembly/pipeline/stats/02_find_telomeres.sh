#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 4 --mem 24gb  --out logs/find_telomeres.log

module load parallel

mkdir -p telomere_reports
ls genomes/*.fasta | parallel -j 4 python  scripts_Hiltunen/find_telomeres.py {} \> telomere_reports/{/.}.telomere_report.txt
