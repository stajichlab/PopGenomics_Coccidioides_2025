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
SBTTEMPLATE=lib/sbt/cocci_pangenome.sbt
IFS=,
tail -n +2 $SAMPLES | sed -n ${N}p | while read RUNACC STRAIN BIOSAMPLE CENTER EXPERIMENT PROJECT ORGANISM FILEBASE NOTES LOCUSTAG
do
    if [[ "$NOTES" == "TooLow" ]]; then
	    echo "skipping $N ($ID) as it is too low coverage ($NOTES)"
	    continue
    elif [[ "$NOTES" == "No RNA" ]]; then
	    echo "skipping $N ($ID) as no good RNA matches ($NOTES)"
	    continue
    fi

    echo "$ID $BASE $SRA $SPECIES $STRAIN"
    SPECIESSTRAINNOSPACE=$(echo -n "$SPECIES $STRAIN" | perl -p -e 's/[\(\)\s]+/_/g')
    SPECIESNOSPACE=$(echo -n "$SPECIES" | perl -p -e 's/[\(\)\s]+/_/g')
    name=$SPECIESSTRAINNOSPACE
    echo "Species is $SPECIESNOSPACE"

    GENOME=$INDIR/${SPECIESSTRAINNOSPACE}.AAFTF.masked.fasta

    if [ -z "$LOCUSTAG" ]; then
        LOCUSTAG=$(echo -n $STRAIN | perl -p -e 's/[\s_\.\-]+//g')
    fi

    funannotate annotate -i $OUTDIR/${name} --cpus $CPUS  \
		--species "$SPECIES" --strain $STRAIN --sbt $SBTTEMPLATE \
	        --busco_db $BUSCODB --rename $LOCUSTAG
		    #$SBTTEMPLATE/$BASE.sbt \
done


