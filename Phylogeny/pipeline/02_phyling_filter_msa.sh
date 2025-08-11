#!/usr/bin/bash -l
#SBATCH -c 24 --mem 96gb --out logs/filter_aln.%A.log
module load phyling
CPU=${SLURM_CPUS_ON_NODE}
if [ -z $CPU ]; then
    CPU=1
fi

marker=ascomycota_odb10
shortmarker=$(basename $marker _odb10)
TYPE=cds
OUTDIR=results
IN=$OUTDIR/align_${TYPE}_${shortmarker}
OUT=$OUTDIR/msa_filter_${TYPE}_${shortmarker}

phyling filter -I $IN -t $CPU -o $OUT --verbose -n 50
