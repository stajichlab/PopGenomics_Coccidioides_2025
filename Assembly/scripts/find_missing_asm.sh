#!/bin/ksh
#SBATCH -p short -c 1 -n 1 --nodes 1 --mem 2G --out logs/find_missing_mask.log

INDIR=asm/AAFTF
OUTDIR=genomes
MASKDIR=RepeatMasker_run
SAMPLES=samples.csv

TORUN=()
IFS=,
N=0
IFS=,
tail -n +2 $SAMPLES | while read RUNACC STRAIN BIOSAMPLE CENTER EXPERIMENT PROJECT ORGANISM FILEBASE NOTES LOCUSTAG
do
    N=$(expr $N + 1)
    for type in AAFTF
    do
    	name=$STRAIN
	OUTNAME=$OUTDIR/${name}.${type}.fasta

	if [[ "$NOTES" == "TooLow" ]]; then
	    #	    echo "skipping $N ($ID $STRAIN) as it is too low coverage ($NOTES)"
	    echo "  skipping $ID $STRAIN too low $NOTES"
	elif [ ! -f $INDIR/${name}.sorted.fasta ]; then
		echo "Cannot find $name.sorted.fasta in $INDIR - may not have been run yet"
	elif [[ -f $OUTNAME && $INDIR/${name}.sorted.fasta -nt $OUTNAME ]]; then
	    echo "need to remove $OUTNAME as $INDIR/${name}.sorted.fasta is newer $N"
	elif [ ! -f $OUTNAME ]; then
	    echo "need to run mask on $N for $OUTNAME"
	    TORUN+=($N)
	fi
    done
done

RUNSET=$(echo "${TORUN[@]}" | perl -p -e 's/ /,/g')

echo "sbatch --array=$RUNSET pipeline/annotation/01_mask.sh"
