#!/bin/bash -l
#SBATCH -p epyc --time 6-0:00:00 -c 24 -n 1 -N 1 --mem 256G --out logs/annotate_train.%a.log

# moved module loading to bottom so we can speed up in case we don't need to actually run this
# module unload miniconda3
# module load funannotate

MEM=256G
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=genomes_to_annotate
ODIR=annotation
SAMPLES=samples.csv
RNAFOLDER=lib/RNASeq
TRAININGCACHE=$(realpath lib/prediction_support/training_cached)
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

export PASAHOME=$HOME/.pasa
echo $PASAHOME
IFS=,
tail -n +2 $SAMPLES | sed -n ${N}p | while read ID BASE SRA SPECIES STRAIN LOCUSTAG BIOPROJECT BIOSAMPLE NOTES
do
    if [[ "$NOTES" == "Too Low" ]]; then
	    echo "skipping $N ($ID) as it is too low coverage ($NOTES)"
	    continue
    fi
    echo "$ID $BASE $SRA $SPECIES $STRAIN"
    SPECIESSTRAINNOSPACE=$(echo -n "$SPECIES $STRAIN" | perl -p -e 's/[\(\)\s]+/_/g')
    SPECIESNOSPACE=$(echo -n "$SPECIES" | perl -p -e 's/[\(\)\s]+/_/g')
    name=$SPECIESSTRAINNOSPACE
    echo "Species is $SPECIESNOSPACE and RNASeq would be $RNAFOLDER/${SPECIESNOSPACE}_R1.fastq.gz"
    
    # previous we were running flye and canu    
    MASKED=$INDIR/${SPECIESSTRAINNOSPACE}.AAFTF.masked.fasta
    echo "in is $MASKED ($INDIR/${name})"
    if [ ! -f $MASKED ]; then
	    echo "no masked file $MASKED"
	    exit
    fi
    if [[ -d $ODIR/${name}/training/genome.fasta && $MASKED -nt $ODIR/${name}/training/genome.fasta ]]; then
	echo "existing training is OLDER than the new genome assembly $MASKED, need to rebuild"
	md5sum $ODIR/${name}/training/genome.fasta
	md5sum $MASKED
	ls -l $MASKED $ODIR/${name}/training/genome.fasta
	exit
    fi
    mkdir -p $ODIR/${name}/training
    if [ -d $TRAININGCACHE/${SPECIESNOSPACE} ]; then
	for nm in trimmomatic normalize
	do
	    if [ ! -e $ODIR/${name}/training/$nm ]; then
		    echo "linking $nm in $ODIR/${name}/training/$nm"
		    ln -s $TRAININGCACHE/${SPECIESNOSPACE}/$nm $ODIR/${name}/training/$nm
	    fi
	done
    fi
    #  load the modules only at bottom for speed since we might skip
    module unload miniconda3
    module load funannotate
    if [[ -f $ODIR/${name}/training/transcript.alignments.gff3 && $MASKED -nt $ODIR/${name}/training/transcript.alignments.gff3  ]]; then
        echo "already generated alignments but  $MASKED is newer than $ODIR/${name}/training/transcript.alignments.gff3, need to remove and rerun"
        exit
    fi
    if [ -f $ODIR/${name}/training/transcript.alignments.gff3 ]; then
        echo "transcript alignments already generated for $name ... skipping"
        exit
    fi
    echo "using $RNAFOLDER/${SPECIESNOSPACE}_R1.fastq.gz and $RNAFOLDER/${SPECIESNOSPACE}_R2.fastq.gz as input RNAseq"
    
    if [[ -f $RNAFOLDER/${SPECIESNOSPACE}_R1.fastq.gz ]]; then
    	funannotate train -i $MASKED -o $ODIR/${name} \
   	     --jaccard_clip --species "$SPECIES" --isolate $STRAIN \
  	     --cpus $CPU --memory ${MEM} \
  	     --left $RNAFOLDER/${SPECIESNOSPACE}_R1.fastq.gz \
	     --right $RNAFOLDER/${SPECIESNOSPACE}_R2.fastq.gz \
  	     --pasa_db mysql
    elif [[ -f $RNAFOLDER/${SPECIESNOSPACE}.fastq.gz ]]; then
	    funannotate train -i $MASKED -o $ODIR/${name} \
   	     --jaccard_clip --species "$SPECIES" --isolate $STRAIN \
  	     --cpus $CPU --memory ${MEM} \
  	     --single $RNAFOLDER/${SPECIESNOSPACE}.fastq.gz \
  	     --pasa_db mysql
    else
	    echo "no RNAfiles for $SPECIESNOSPACE in $RNAFOLDER"
	    exit
    fi
done
