#!/usr/bin/bash -l
#SBATCH -p short -c 24 --mem 16gb --out logs/map_stats.log

module load samtools
CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
pushd aln
parallel -j $CPU samtools flagstat {} \> {.}.flagstat ::: $(ls *.cram)
grep mapped *.flagstat | grep -v -P 'mate|primary' | perl -p -e 's/\.flagstat:/\t/; s/\s+\+\s+\d+\s+mapped\s+\((\d+\.\d+)\%\s+.+/\t$1/' | sort -k3,3nr > ../aln.mapped.tsv
