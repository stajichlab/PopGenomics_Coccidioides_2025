#!/usr/bin/bash -l
#SBATCH -p short -c 8 --mem 24gb --out logs/organize_by_enzymeclass.log

CPU=1
if [ ! -z $SLURM_CPUS_PER_TASK ]; then
    CPU=$SLURM_CPUS_PER_TASK
fi
INFOLDER=results/bigscape_v2/output_files/protocl_simplematch_hybridoff_category_2025-08-10_18-39-53_c0.3
OUTFOLDER=results/bigscape_v2.seq_by_class
INFASTA=results/antismash_regions/fasta
mkdir -p $OUTFOLDER

python scripts/organize_clusterseq.py \
-a $INFOLDER/record_annotations.tsv --threads $CPU \
$INFASTA $OUTFOLDER
