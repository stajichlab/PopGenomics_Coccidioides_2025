#!/usr/bin/bash -l
#SBATCH -c 24 --mem 96gb -p short --out logs/make_tree_cds.%A.log
module load phyling
CPU=${SLURM_CPUS_ON_NODE}
if [ -z $CPU ]; then
    CPU=1
fi

marker=ascomycota_odb10
shortmarker=$(basename $marker _odb10)
TYPE=cds
OUTDIR=results
IN=$OUTDIR/msa_filter_${TYPE}_${shortmarker}
OUT=$OUTDIR/tree_${TYPE}_${shortmarker}

phyling tree -I $IN -M ft -t $CPU -o $OUT --verbose
