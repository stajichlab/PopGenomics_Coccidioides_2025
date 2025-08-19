#!/usr/bin/bash -l
#SBATCH --time 3-0:00:00 -c 16 -N 1 -n 1 --mem 16G --out logs/annotate_predict_noCQ.%a.log

module load funannotate
which augustus
# this will define $SCRATCH variable if you don't have this on your system you can basically do this depending on
# where you have temp storage space and fast disks
module load workspace/scratch

CPU=1
if [ ! -z $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

BUSCO=onygenales_odb10 # This could be changed to the core BUSCO set you want to use
INDIR=scaffolded_genomes
OUTDIR=annotation
mkdir -p $OUTDIR
SAMPLES=samples.csv

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
MAX=$(wc -l $SAMPLES | awk '{print $1}')

if [ $N -gt $MAX ]; then
    echo "$N is too big, only $MAX lines in $SAMPLES"
    exit
fi

export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/3.5/config)
export FUNANNOTATE_DB=/bigdata/stajichlab/shared/lib/funannotate_db

SEED_SPECIES=coccidioides_immitis
SEQCENTER=UCR
IFS=,
tail -n +2 $SAMPLES | sed -n ${N}p | while read RUNACC STRAIN BIOSAMPLE CENTER EXPERIMENT PROJECT ORGANISM FILEBASE NOTES LOCUSTAG
do
    if [[ "$NOTES" == "TooLow" ]]; then
	echo "skipping $N ($ID) as it is too low coverage ($NOTES)"
	continue
    fi
    if [[ "$NOTES" == "Skip" ]]; then
	    echo "Skipping $N ($ID) as it was marked to be skipped"
	    continue
    fi
    SPECIESSTRAINNOSPACE=$(echo -n "$ORGANISM" | perl -p -e 's/[\(\)\s]+/_/g')
    SPECIESNOSPACE=$(echo -n "$ORGANISM" | perl -p -e 's/[\(\)\s]+/_/g')
    name=$STRAIN
    MASKED=$INDIR/${name}.masked.fasta
    echo "masked is $MASKED ($INDIR/${name}.masked.fasta)"
    if [ ! -f $MASKED ]; then
        echo "no masked file $MASKED"
           exit
    fi
    if [ -z "$LOCUSTAG" ]; then
	    LOCUSTAG=$(echo -n $STRAIN | perl -p -e 's/[\s_\.\-]+//g')
    fi
    echo "LOCUS is $LOCUSTAG MASKED is $MASKED"
    if [[ -f $OUTDIR/${name}/predict_results/${SPECIESNOSPACE}_$STRAIN.proteins.fa && $OUTDIR/${name}/predict_results/${SPECIESNOSPACE}_${STRAIN}.proteins.fa -nt $MASKED ]]; then
	    echo "Already run prediction step, skipping, to force remove the proteins.fa file"
    fi
    if [[ -f $OUTDIR/${name}/predict_misc/protein_alignments.gff3 && $MASKED -nt $OUTDIR/${name}/predict_misc/protein_alignments.gff3 ]]; then
	    echo "$MASKED is newer than $OUTDIR/${name}/predict_misc/protein_alignments.gff3, need to remove existing files to ensure clean re-run"
	    exit
    fi
    # we cannot re-rerun if this file is in place
    rm -rf  $OUTDIR/${name}/predict_misc/CodingQuarry
    if [[ -f $OUTDIR/${name}/predict_misc/augustus.gff3 && $OUTDIR/${name}/training/funannotate_train.stringtie.gtf -nt $OUTDIR/${name}/predict_misc/augustus.gff3 ]]; then
	    rm -rf $OUTDIR/${name}/predict_misc/augustus.* $OUTDIR/${name}/predict_misc/snap* $OUTDIR/${name}/predict_misc/EVM \
		    $OUTDIR/${name}/predict_misc/hints* $OUTDIR/${name}/predict_misc/genemark.* $OUTDIR/${name}/predict_misc/gene_predictions.gff3
    fi
    time funannotate predict --cpus $CPU --keep_no_stops --SeqCenter $SEQCENTER \
		--busco_db $BUSCO --header_length 24 \
		--strain $STRAIN --min_training_models 40 \
		--AUGUSTUS_CONFIG_PATH $AUGUSTUS_CONFIG_PATH \
		-i $MASKED --name $LOCUSTAG --max_intronlen 1500 \
		--protein_evidence $FUNANNOTATE_DB/uniprot_sprot.fasta \
		-s "$ORGANISM" -o $OUTDIR/${name}  --tmpdir $SCRATCH \
		--busco_seed_species $SEED_SPECIES -w codingquarry:0
done
