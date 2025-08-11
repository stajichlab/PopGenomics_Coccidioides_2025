#!/usr/bin/bash -l
#SBATCH -c 96 -N 1 -n 1 --mem 128gb --out logs/phyling_align.log

module load phyling
CPU=96
if  [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi
marker=ascomycota_odb10
shortmarker=$(basename $marker _odb10)
IN=input
TYPE=cds
OUTDIR=results
OUT=${OUTDIR}/align_${TYPE}_${shortmarker}
mkdir -p $OUT
phyling align -I $IN -m $marker -t $CPU -o $OUT -v
