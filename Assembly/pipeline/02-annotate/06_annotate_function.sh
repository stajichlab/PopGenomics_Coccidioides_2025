#!/bin/bash -l
#SBATCH --nodes 1 -c 24 -n 1 --mem 64G --out logs/annotate_function.%a.log
# note this doesn't need that much memory EXCEPT for the XML -> tsv parsing that happens when you provided an interpro XML file

MEM=64G
module load funannotate

CPUS=$SLURM_CPUS_ON_NODE

if [ ! $CPUS ]; then
    CPUS=2
fi
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

INDIR=scaffolded_genomes
OUTDIR=annotation
BUSCODB=onygenales_odb10
SBTTEMPLATE=lib/sbt/Cocci.sbt
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
    # need to skip Contam too

    echo "Species is $ORGANISM"

    if [ -z "$LOCUSTAG" ]; then
	    echo "no LOCUSTAG for $STRAIN"
        LOCUSTAG=$(echo -n $STRAIN | perl -p -e 's/[\s_\.\-]+//g')
    fi
    name=$STRAIN
    ANTISMASH=$OUTDIR/$name/antismash_local/$SPECIESNOSPACE.gbk
    ARGS=()
    if [[ -d $(dirname $ANTISMASH) && -s $ANTISMASH ]]; then
	ARGS+=(--antismash $ANTISMASH)
    fi

    funannotate annotate -i $OUTDIR/${name} --cpus $CPUS --header_length 24 \
		--species "$ORGANISM" --strain $STRAIN --sbt $SBTTEMPLATE \
	        --busco_db $BUSCODB --rename $LOCUSTAG "${ARGS[@]}" 
		   #--rename $LOCUSTAG
done
