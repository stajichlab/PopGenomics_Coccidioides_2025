#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 24 --mem 64gb --out logs/genomescope.%a.log -a 1-186

module load workspace/scratch
module load samtools
module load jellyfish
module load R

GENOMESCOPE=genomescope

CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
  CPU=$SLURM_CPUS_ON_NODE
fi
N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
  N=$1
fi
if [ -z $N ]; then
  echo "cannot run without a number provided either cmdline or --array in sbatch"
  exit
fi
FASTQFOLDER=input
SAMPLEFILE=samples.csv
MAX=$(wc -l $SAMPLEFILE | awk '{print $1}')
if [ $N -gt $MAX ]; then
  echo "$N is too big, only $MAX lines in samplefile=$SAMPLEFILE"
  exit
fi
mkdir -p $GENOMESCOPE
JELLYFISHSIZE=1000000000
IFS=,
KMER=21
READLEN=150 # note this assumes all projects are 150bp reads which they may not be
IFS=, # set the delimiter to be ,
tail -n +2 $SAMPLEFILE | sed -n ${N}p | while read ID BASE SRA SPECIES STRAIN LOCUSTAG BIOPROJECT BIOSAMPLE NOTES LOCUSTAG
do
    # for this project two different file base but this needs fixing otherse)
  if [ ! -s $GENOMESCOPE/$STRAIN.histo ]; then
    echo "running genomescope on $STRAIN"
	  if [ ! -f $FASTQFOLDER/${BASE}_1.fastq.gz ]; then
      COUNT=$(echo -n $BASE | perl -p -e 's/;/\n/g' | wc -l | awk '{print $1}')
      if [[ $COUNT == 0 ]]; then
        FNAME=$FASTQFOLDER/${BASE}
        echo "FNAME is $FNAME"
      else
        for file in $(echo -n $BASE | perl -p -e 's/;/,/g')
        do
          cat $FASTQFOLDER/$file
        done > $SCRATCH/${ID}_all.fastq.gz
        FNAME=$SCRATCH/${ID}_all.fastq.gz
      fi
      jellyfish count -C -m $KMER -s $JELLYFISHSIZE -t $CPU -o $SCRATCH/$STRAIN.jf <(pigz -dc $FNAME)
      rm -f $SCRATCH/${ID}_all.fastq.gz
    else
      jellyfish count -C -m $KMER -s $JELLYFISHSIZE -t $CPU -o $SCRATCH/$STRAIN.jf <(pigz -dc $FASTQFOLDER/${BASE}_[12].fastq.gz)
	  fi
    jellyfish histo -t $CPU $SCRATCH/$STRAIN.jf > $GENOMESCOPE/$STRAIN.histo
    Rscript scripts/genomescope.R $GENOMESCOPE/$STRAIN.histo $KMER $READLEN $GENOMESCOPE/$STRAIN/
  else
    echo "skipping ... $STRAIN already processed"
  fi
done
