#!/bin/ksh
#SBATCH -p short logs/find_missing_masked.log

CPU=1

INDIR=scaffolded_genomes
SAMPLES=samples.csv
IFS=,
N=1
mkdir -p empty
TORUN=()
IFS=,
tail -n +2 $SAMPLES | while read RUNACC STRAIN BIOSAMPLE CENTER EXPERIMENT PROJECT ORGANISM FILEBASE NOTES LOCUSTAG
do
    if [[ "$NOTES" == "TooLow" || "$NOTES" == "Contaminant" ]]; then
        echo "skipping $N ($ID) as low coverage or contaminant ($NOTES)" >> check_trained_msg.log
        N=$(expr $N + 1)
        continue
    fi
    SPECIESSTRAINNOSPACE=$(echo -n "$ORGANISM $STRAIN" | perl -p -e 's/[\(\)\s]+/_/g')
    SPECIESNOSPACE=$(echo -n "$ORGANISM" | perl -p -e 's/[\(\)\s]+/_/g')
    name=$STRAIN
    proteins=annotation/${name}/predict_results/${SPECIESSTRAINNOSPACE}.proteins.fa
    if [ ! -f $INDIR/$name.masked.fasta ]; then
	echo -e "\tCannot find $INDIR/${name}.masked.fasta in $INDIR - may not have been run yet ($N)" 1>&2
    elif [ ! -f $proteins ]; then
        echo "need to run annotate on $name ($N) no $proteins" 1>&2
	if [ ! -s annotation/${name}/predict_results/augustus.gff3 ]; then
	    echo "echo annotation/${name}" >>  delete_$$.sh
	    echo "/usr/bin/rm -rf annotation/${name}/predict_misc/busco*" >> delete_$$.sh
	    echo "mv annotation/${name}/predict_misc/EVM_busco annotation/${name}/predict_misc/EVM_busco.b" >> delete_$$.sh
	    echo "/usr/bin/rm -rf annotation/${name}/predict_misc/hints.*" >> delete_$$.sh
	fi
	if [ ! -s annotation/${name}/predict_misc/genemark/genemark.gtf ]; then
		 echo "rm -rf annotation/${name}/predict_misc/genemark*" >> delete_$$.sh
	fi
	TORUN+=($N)
    fi
    N=$(expr $N + 1)
done

echo 'for file in annotate/*/predict_misc/EVM_busco.b; do rsync -a --delete ./empty/ $file/; rmdir $file; done' >> delete_$$.sh
RUNSET=$(echo "${TORUN[@]}" | perl -p -e 's/ /,/g')

echo "sbatch --array=$RUNSET pipeline/annotation/03_predict.sh"
