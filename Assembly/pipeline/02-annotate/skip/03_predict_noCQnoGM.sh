#!/usr/bin/bash -l
#SBATCH --time 3-0:00:00 --ntasks 16 --nodes 1 --mem 24G --out logs/annotate_predict.%a.log

module load funannotate

# this will define $SCRATCH variable if you don't have this on your system you can basically do this depending on
# where you have temp storage space and fast disks
module load workspace/scratch

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

BUSCO=basidiomycota_odb10 # This could be changed to the core BUSCO set you want to use
INDIR=genomes_to_annotate
OUTDIR=annotation
mkdir -p $OUTDIR
SAMPFILE=samples.csv

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "need to provide a number by --array or cmdline"
        exit
    fi
fi
MAX=$(wc -l $SAMPFILE | awk '{print $1}')

if [ $N -gt $MAX ]; then
    echo "$N is too big, only $MAX lines in $SAMPFILE"
    exit
fi

#export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/3.3/config)
export FUNANNOTATE_DB=/bigdata/stajichlab/shared/lib/funannotate_db

SEED_SPECIES=ustilago
SEQCENTER=UCR
IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read ID BASE SRA SPECIES STRAIN LOCUSTAG BIOPROJECT BIOSAMPLE NOTES
do
    if [[ "$NOTES" == "Too Low" ]]; then
	echo "skipping $N ($ID) as it is too low coverage ($NOTES)"
	continue
    fi
    SPECIESSTRAINNOSPACE=$(echo -n "$SPECIES $STRAIN" | perl -p -e 's/[\(\)\s]+/_/g')
    SPECIESNOSPACE=$(echo -n "$SPECIES" | perl -p -e 's/[\(\)\s]+/_/g')
    name=$SPECIESSTRAINNOSPACE
    MASKED=$INDIR/${name}.AAFTF.masked.fasta
    echo "masked is $MASKED ($INDIR/${name}.masked.fasta)"
    if [ ! -f $MASKED ]; then
        echo "no masked file $MASKED"
           exit
    fi
    if [ -z "$LOCUSTAG" ]; then
	    LOCUSTAG=$(echo -n $STRAIN | perl -p -e 's/[\s_\.\-]+//g')
    fi
    echo "LOCUS is $LOCUSTAG MASKED is $MASKED"
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
		--busco_db $BUSCO --optimize_augustus -w codingquarry:0 genemark:0 \
		--strain $STRAIN --min_training_models 100 \
		--AUGUSTUS_CONFIG_PATH $AUGUSTUS_CONFIG_PATH \
		-i $MASKED --name $LOCUSTAG \
		--protein_evidence $FUNANNOTATE_DB/uniprot_sprot.fasta \
		-s "$SPECIES" -o $OUTDIR/${name}  --tmpdir $SCRATCH
		#--busco_seed_species $SEED_SPECIES
done
