#!/usr/bin/bash -l

#SBATCH -N 1 -n 1 -c 24 --mem 24G -J readcount --out logs/bbcount.%a.log --time 48:00:00
module load BBMap
module load workspace/scratch

hostname
MEM=24
CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ ! $N ]; then
    N=$1
    if [ ! $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

INDIR=input
FASTQDIR=$INDIR
SAMPFILE=samples.csv
TEMP=$SCRATCH
ASM=genomes
OUTDIR=$(realpath mapping_report)
mkdir -p $OUTDIR

IFS=, # set the delimiter to be ,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read RUNACC STRAIN BIOSAMPLE CENTER EXPERIMENT PROJECT ORGANISM FILEBASE NOTES LOCUSTAG
do
    if [[ "$NOTES" == "TooLow" ]]; then
		echo "skipping $N ($ID/$STRAIN) as it is too low coverage ($NOTES)"
    fi
    DONETHIS=0
    for type in AAFTF
    do
		REPORTOUT=${STRAIN}.${type}
		SORTED=$(realpath $ASM/${STRAIN}.${type}.fasta)
		if [[ -s $OUTDIR/${REPORTOUT}.bbmap_summary.txt && $OUTDIR/${REPORTOUT}.bbmap_summary.txt -nt $SORTED ]]; then
			DONETHIS=1
		else
			DONETHIS=0
			break
		fi
    done
    
    if [[ $DONETHIS == "1" ]]; then
		exit
    fi
	
    LEFTIN=$(ls $FASTQDIR/${RUNACC}${FILEBASE} | sed -n 1p)
    RIGHTIN=$(ls $FASTQDIR/${RUNACC}${FILEBASE} | sed -n 2p)
    LEFT=$(realpath $LEFTIN)
    RIGHT=$(realpath $RIGHTIN)
    echo "$LEFT $RIGHT"
    for type in AAFTF
    do
		SORTED=$(realpath $ASM/${STRAIN}.${type}.fasta)
		REPORTOUT=${STRAIN}.${type}
		if [ -s $SORTED ]; then
			pushd $SCRATCH
			if [ ! -s $OUTDIR/${REPORTOUT}.bbmap_covstats.txt ]; then
				bbmap.sh -Xmx${MEM}g ref=$SORTED in=$LEFT in2=$RIGHT covstats=$OUTDIR/${REPORTOUT}.bbmap_covstats.txt  statsfile=$OUTDIR/${REPORTOUT}.bbmap_summary.txt
			fi
			popd
		else
			echo "no $SORTED for $STRAIN ($ID) file"
		fi
    done
done

