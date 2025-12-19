#!/bin/bash -l
#SBATCH -p epyc -c 8 -n 1 --nodes 1 --mem 48G --out logs/annotate_noscaff_mask.%a.log

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=asm/AAFTF
OUTDIR=unscaffolded_genomes
MASKDIR=RepeatMasker_run_noscaffold
SAMPLES=samples.csv
RMLIBFOLDER=lib/repeat_library
mkdir -p $RMLIBFOLDER $MASKDIR $OUTDIR
RMLIBFOLDER=$(realpath $RMLIBFOLDER)

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

IFS=,
tail -n +2 $SAMPLES | sed -n ${N}p | while read RUNACC STRAIN BIOSAMPLE CENTER EXPERIMENT PROJECT ORGANISM FILEBASE NOTES LOCUSTAG
do
    if [[ "$NOTES" == "TooLow" ]]; then
		echo "skipping $N ($ID $STRAIN) as it is too low coverage ($NOTES)"
		continue
    fi
    SPECIESNOSPACE=$(echo -n "${ORGANISM}" | perl -p -e 's/[\(\)\s]+/_/g')
    name=$STRAIN
    ASMFILE=$INDIR/$STRAIN.sorted.fasta
    if [ ! -f $ASMFILE ]; then
		echo "Cannot find $ASMFILE in $INDIR - may not have been run yet"
		exit
    fi
    OUTNAME=$OUTDIR/$STRAIN.masked.fasta
    if [[ ! -s $OUTNAME || $ASMFILE -nt $OUTNAME ]]; then
	mkdir -p $MASKDIR/${name}
	GENOME=$(realpath $ASMFILE)
	LIBRARY=$RMLIBFOLDER/$SPECIESNOSPACE.repeatmodeler.lib
	if [[ ! -f $LIBRARY ]]; then
		echo "cannot find the RepeatModeler library for $SPECIESNOSPACE"
		exit
#	    module load RepeatModeler
#	    pushd $MASKDIR/${name}
#	    rm -rf RM_*
#	    BuildDatabase -name $STRAIN $GENOME
#	    RepeatModeler -threads $CPU -database $STRAIN -LTRStruct
#	    rsync -a RM_*/consensi.fa.classified $LIBRARY
#	    rsync -a RM_*/families-classified.stk $RMLIBFOLDER/$SPECIESNOSPACE.repeatmodeler.stk
#	    popd
	fi
	if [ -f $LIBRARY ]; then
		module load RepeatMasker
		RepeatMasker -e ncbi -xsmall -s -pa $CPU -lib $LIBRARY -dir $MASKDIR/${name} -gff $ASMFILE 
	fi	
	rsync -a $MASKDIR/${name}/$STRAIN.sorted.fasta.masked $OUTNAME
    else
	echo "Skipping ${name} as masked file already exists and is newer than current assembly file"
    fi
done
