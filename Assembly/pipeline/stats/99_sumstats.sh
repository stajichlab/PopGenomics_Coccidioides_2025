#!/usr/bin/bash -l
#SBATCH -p short

perl scripts/asm_stats_scaffolds.pl scaffolded_genomes > asm_scaffolded_stats.tsv
perl scripts/asm_stats.pl genomes > asm_stats.tsv
