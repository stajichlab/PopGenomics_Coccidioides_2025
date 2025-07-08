#!/usr/bin/bash -l
#SBATCH -p short

echo -e "STRAIN\tMIN_HET\tMAX_HET" > heterozygocity.tsv
pushd genomescope
grep Het */summary.txt | perl -p -e 's/\/summary.txt:\S+\s+/\t/s; s/\%//g; s/ +/\t/g; s/\s+$/\n/g;' >> ../heterozygocity.tsv

