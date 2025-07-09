#!/usr/bin/bash
#SBATCH -n 1 -c 64 -N 1 -p short --mem 72G --out logs/mosdepth.parallel.log
#SBATCH -J modepth
# This one goes to 11
module load parallel
CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
 CPU=2
fi
module load mosdepth
mkdir -p coverage/mosdepth
source config.txt

for WINDOW in 5000 10000 50000
do
	parallel --jobs $CPU mosdepth -f $REFGENOME -T 1,10,50,100,200 -n --by $WINDOW -t 2 "{= s:$ALNFOLDER\/:coverage/mosdepth/:; s:\.$HTCEXT:.${WINDOW}bp: =}" {} ::: $ALNFOLDER/*.$HTCEXT
done

mkdir -p plots
Rscript scripts/plot_mosdepth_CNV.R

