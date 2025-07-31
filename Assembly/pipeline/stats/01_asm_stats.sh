#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 2 --mem 4gb --out logs/stats.log

SAMPFILE=samples.csv
INDIR=asm
OUTDIR=genomes

mkdir -p $OUTDIR
module load AAFTF

IFS=, # set the delimiter to be ,
tail -n +2 $SAMPFILE | while read RUNACC STRAIN BIOSAMPLE CENTER EXPERIMENT PROJECT ORGANISM FILEBASE NOTES LOCUSTAG
do
	ID=$STRAIN
    if [[ "$NOTES" == "TooLow" ]]; then
		echo "skipping $N ($ID) as it is too low coverage ($NOTES)"
		continue
    fi
    # this used to be setup for different assembly methods too
    for type in AAFTF
    do
	if [ ! -f $INDIR/$type/$ID.sorted.fasta ]; then
		echo "cannot find $INDIR/$type/$ID.sorted.fasta"
		continue
	fi
	if [[ ! -f $OUTDIR/$STRAIN.$type.fasta || $INDIR/$type/$ID.sorted.fasta -nt $OUTDIR/$STRAIN.$type.fasta ]]; then
		echo "need to re-copy or copy 1st time $INDIR/$type/$ID.sorted.fasta $OUTDIR/$STRAIN.$type.fasta"
	fi
	rsync -avL $INDIR/$type/$ID.sorted.fasta $OUTDIR/$STRAIN.$type.fasta
	if [[ ! -f $OUTDIR/$STRAIN.$type.stats.txt || $OUTDIR/$STRAIN.$type.fasta -nt $OUTDIR/$STRAIN.$type.stats.txt ]]; then
		AAFTF assess -i $OUTDIR/$STRAIN.$type.fasta -r $OUTDIR/$STRAIN.$type.stats.txt
	fi
    done
done
