#!/usr/bin/bash -l
#SBATCH -p short

perl scripts/asm_stats_scaffolds.pl scaffolded_genomes > asm_scaffolded_stats.$$.tsv

(head -n 1 asm_scaffolded_stats.$$.tsv  && tail -n +2 asm_scaffolded_stats.$$.tsv | sort -k3,3nr) > asm_scaffolded_stats.tsv

perl scripts/asm_stats.pl genomes > asm_stats.$$.tsv
(head -n 1 asm_stats.$$.tsv  && tail -n +2 asm_stats.$$.tsv | sort -k3,3nr) > asm_stats.tsv
rm -f asm_stats.$$.tsv asm_scaffolded_stats.$$.tsv
