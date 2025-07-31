#!/bin/bash -l
#SBATCH -p short -c 1 -n 1 -N 1 --mem 1G 

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

INDIR=scaffolded_genomes
ODIR=annotation
SAMPLES=samples.csv
N=1
rm -f check_trained_msg.log to_run_training.txt
IFS=,
tail -n +2 $SAMPLES | while read RUNACC STRAIN BIOSAMPLE CENTER EXPERIMENT PROJECT ORGANISM FILEBASE NOTES LOCUSTAG
do
    if [[ "$NOTES" == "TooLow" ]]; then
        echo "skipping $N ($ID) as it is too low coverage ($NOTES)" >> check_trained_msg.log
        N=$(expr $N + 1)
        continue
    fi
    SPECIESNOSPACE=$(echo -n "$ORGANISM" | perl -p -e 's/[\(\)\s]+/_/g')
    name=$STRAIN

    MASKED=$INDIR/${name}.masked.fasta

    if [ ! -f $MASKED ]; then
        echo "no masked file $MASKED"
        N=$(expr $N + 1)
        continue
    fi
    if [[ -f $ODIR/${name}/training/funannotate_train.pasa.gff3 && $MASKED -nt $ODIR/${name}/training/funannotate_train.pasa.gff3 ]]; then
        echo "already generated alignments but  $MASKED is newer than $ODIR/${name}/training/funannotate_train.pasa.gff3, need to remove and rerun"
        N=$(expr $N + 1)
    fi
    if [ -f $ODIR/${name}/training/funannotate_train.pasa.gff3 ]; then
        echo "transcript alignments already generated for $name ($ODIR/${name}/training/trinity.alignments.gff3) ... skipping" >> check_trained_msg.log
        N=$(expr $N + 1)
        continue
    fi
    echo "$N $ODIR/${name} needs to be run" 
    echo -n "$N," >> to_run_training.txt
    N=$(expr $N + 1)
done
