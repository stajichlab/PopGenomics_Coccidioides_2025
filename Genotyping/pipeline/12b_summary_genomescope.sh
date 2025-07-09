#!/usr/bin/bash -l
#SBATCH -p short
mkdir -p summary_stats
OUTFILE=summary_stats/heterozygocity.tsv
echo -e "STRAIN\tMIN_HET\tMAX_HET" > $OUTFILE
pushd genomescope
grep Het */summary.txt | perl -p -e 's/\/summary.txt:\S+\s+/\t/s; s/\%//g; s/ +/\t/g; s/\s+$/\n/g;' >> ../$OUTFILE

