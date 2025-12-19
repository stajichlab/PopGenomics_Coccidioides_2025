#!/bin/bash -l
#SBATCH -p epyc --time 5-0:00:00 -c 16 -n 1 -N 1 --mem 24G --out logs/annotate_noscaffold_train.%a.log

module load singularity
HOSTNAME=$(hostname -s)
echo "Running on $HOSTNAME"
# Define stop mysqldb
RUNID=$$
stop_mysqldb() { singularity instance stop mysqldb$RUNID; }

# Define error handler
error_exit()
{
    stop_mysqldb
	echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2
	exit 1
}

# Set trap to ensure mysqldb is stopped
trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM

# Set some vars
mkdir -p $SCRATCH/db $SCRATCH/conf
rsync -a ~/bigdata/mysql/db/mysql $SCRATCH/db/ || error_exit "Failed to copy mysql data"
cp  ~/.pasa/pasa_conf/my.cnf $SCRATCH/conf/my.cnf || error_exit "Failed to copy pasa config file"
cp ~/.pasa/pasa_conf/conf.txt $SCRATCH/conf/pasa-local-${HOSTNAME}.config.txt 
PORT=$(shuf -i3000-4999 -n1)

export SINGULARITY_BINDPATH=$SCRATCH
export PASACONF=$SCRATCH/conf/pasa-local-${HOSTNAME}.config.txt
#export PASACONF=$SINGULARITYENV_PASACONF
sed -i "s/^MYSQLSERVER.*$/MYSQLSERVER=${HOSTNAME}:${PORT}/" ${PASACONF}
perl -i -p -e "s/port = \d+/port = ${PORT}/" $SCRATCH/conf/my.cnf
SIF=/bigdata/stajichlab/shared/lib/mariadb/mariadb.sif
# Start Database
singularity instance start --writable-tmpfs -B $SCRATCH/conf/my.cnf:/etc/mysql/my.cnf,$SCRATCH/db/:/var/lib/mysql,$SCRATCH/conf:/usr/conf $SIF mysqldb$RUNID /usr/bin/mysqld_safe
MEM=24G
CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=unscaffolded_genomes
ODIR=unscaffolded_annotation
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

echo $PASAHOME
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

    echo "$STRAIN $ORGANISM"
    SPECIESSTRAINNOSPACE=$(echo -n "$ORGANISM" | perl -p -e 's/[\(\)\s]+/_/g')
    SPECIESNOSPACE=$(echo -n "$ORGANISM" | perl -p -e 's/[\(\)\s]+/_/g')
    name=$STRAIN
    echo "Species is $SPECIESNOSPACE and RNASeq would be $RNAFOLDER/${SPECIESNOSPACE}_R1.fq.gz"
    
    MASKED=$INDIR/${name}.masked.fasta
    echo "input genome is $MASKED ($INDIR/${name})"
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
    module load funannotate
    module load trinity-rnaseq
    export PASAHOME=$HOME/.pasa
    if [[ -f $ODIR/${name}/training/funannotate_train.pasa.gff3 && $MASKED -nt $ODIR/${name}/training/funannotate_train.pasa.gff3 ]]; then
        echo "already generated alignments but  $MASKED is newer than $ODIR/${name}/training/funannotate_train.pasa.gff3, need to remove and rerun"
        exit
    fi
    if [[ -f $ODIR/${name}/training/funannotate_train.pasa.gff3 ]]; then
        echo "transcript alignments already generated for $name ($ODIR/${name}/training/trinity.alignments.gff3) ... skipping"
	continue
    fi
    echo "using $RNAFOLDER/${SPECIESNOSPACE}_R1.fq.gz and $RNAFOLDER/${SPECIESNOSPACE}_R2.fq.gz as input RNAseq"
    
    if [[ -f $RNAFOLDER/${SPECIESNOSPACE}_R1.fq.gz ]]; then
        funannotate train -i $MASKED -o $ODIR/${name} \
   	        --jaccard_clip --species "$ORGANISM" --strain $STRAIN \
  	        --cpus $CPU --memory ${MEM} --header_length 24 \
  	        --left $RNAFOLDER/${SPECIESNOSPACE}_R1.fq.gz \
	        --right $RNAFOLDER/${SPECIESNOSPACE}_R2.fq.gz \
  	        --pasa_db mysql --no-progress --min_coverage 4 
    elif [[ -f $RNAFOLDER/${SPECIESNOSPACE}.fq.gz ]]; then
	    funannotate train -i $MASKED -o $ODIR/${name} \
   	        --jaccard_clip --species "$SPECIES" --strain $STRAIN \
  	        --cpus $CPU --memory ${MEM} --header_length 24 \
  	        --single $RNAFOLDER/${SPECIESNOSPACE}.fq.gz \
  	        --pasa_db mysql --no-progress --min_coverage 4
    else
	    echo "no RNAfiles for $SPECIESNOSPACE in $RNAFOLDER"
	    exit
    fi
done
# shutdown the DB
#stop_mysqldb
exit 0
