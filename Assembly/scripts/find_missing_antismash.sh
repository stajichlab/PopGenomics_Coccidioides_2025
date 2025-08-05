#!/bin/ksh
#SBATCH -p short logs/find_missing_masked.log

CPU=1

INDIR=annotation
SAMPLES=samples.csv
IFS=,
N=1
mkdir -p empty
TORUN=()
IFS=,
tail -n +2 $SAMPLES | while read RUNACC STRAIN BIOSAMPLE CENTER EXPERIMENT PROJECT ORGANISM FILEBASE NOTES LOCUSTAG
do
    if [[ "$NOTES" == "TooLow" ]]; then
        echo "skipping $N ($ID) as it is too low coverage ($NOTES)" >> check_trained_msg.log
        N=$(expr $N + 1)
        continue
    fi
    name=$STRAIN
    if [[ ! -f $INDIR/${name}/antismash_local/index.html ]]; then
	    TORUN+=($N)
    fi
    N=$(expr $N + 1)
done
RUNSET=$(echo "${TORUN[@]}" | perl -p -e 's/ /,/g' | perl -p -e 's/,\s*$//')
echo "sbatch --array=$RUNSET pipeline/annotate/05a_antismash_local.sh"
