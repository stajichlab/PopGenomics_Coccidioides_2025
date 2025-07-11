#!/usr/bin/bash -l
#SBATCH -N 1 -c 16 -n 1 --mem 32gb --out logs/bwa.%a.log --time 2:00:00 -a 1-90 -p short
module load bwa-mem2
module load samtools
module load picard
module load gatk/4.6.1.0
module load java
module load workspace/scratch
MEM=32g

TEMP=$SCRATCH

if [ -f config.txt ]; then
  source config.txt
fi
TEMP=$SCRATCH
if [ -z $REFGENOME ]; then
  echo "NEED A REFGENOME - set in config.txt and make sure 00_index.sh is run"
  exit
fi

if [ ! -f $REFGENOME.dict ]; then
  echo "NEED a $REFGENOME.dict - make sure 00_index.sh is run"
fi
mkdir -p $TMPOUTDIR $ALNFOLDER $UNMAPPED

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

MAX=$(wc -l $SAMPFILE | awk '{print $1}')
if [ $N -gt $MAX ]; then
  echo "$N is too big, only $MAX lines in $SAMPFILE"
  exit
fi

IFS=,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read RUNACC STRAIN BIOSAMPLE CENTER EXPERIMENT PROJECT ORGANISM FILEBASE
do
  FINALFILE=$ALNFOLDER/$STRAIN.$HTCEXT
  echo "To process $STRAIN and produce $FINALFILE"
  if [ ! -s $FINALFILE ]; then
    BAMSTOMERGE=()
    # support multiple FASTQ per strain would be split by ';'
    for BASE in $(echo $RUNACC | perl -p -e 's/\;/,/g');
    do
      # END THIS PART IS PROBABLY PROJECT SPECIFIC
      echo "STRAIN is $STRAIN BASE is $BASE BASE is $BASE"
      TMPBAMFILE=$TEMP/$BASE.unsrt.bam
      SRTED=$TEMP/$BASE.srt.bam
      DDFILE=$TEMP/$BASE.DD.bam

      FINALFILE=$ALNFOLDER/$STRAIN.$HTCEXT
      READGROUP="@RG\tID:$BASE\tSM:$STRAIN\tLB:$BASE\tPL:illumina\tCN:$RGCENTER"
      echo "$TMPBAMFILE $READGROUP"

      if [ ! -s $DDFILE ]; then
        if [ ! -s $SRTED ]; then
          if [ -e $PAIR1 ]; then
            if [ ! -f $SRTED ]; then
              # potential switch this to bwa-mem2 for extra speed
              bwa-mem2 mem -t $CPU -R $READGROUP $REFGENOME $FASTQFOLDER/${BASE}${FILEBASE} > $SCRATCH/${BASE}.sam
              samtools fixmate --threads $CPU -u -O BAM $SCRATCH/${BASE}.sam $SCRATCH/fixmate.bam
              samtools sort --threads $CPU -O BAM -o $SRTED -T $TEMP $SCRATCH/fixmate.bam
            fi
          else
            echo "Cannot find $BASE, skipping $STRAIN"
            exit
          fi
        fi # SRTED file exists or was created by this block

        time java -jar $PICARD MarkDuplicates -I $SRTED -O $DDFILE \
          -METRICS_FILE logs/$STRAIN.dedup.metrics -CREATE_INDEX true -VALIDATION_STRINGENCY SILENT
        if [ -f $DDFILE ]; then
          rm -f $SRTED
        fi
      fi # DDFILE is created after this or already exists
      BAMSTOMERGE+=( $DDFILE )
    done
    samtools merge --write-index -O $HTCFORMAT --threads $CPU --reference $REFGENOME $FINALFILE "${BAMSTOMERGE[@]}"
    if [ -f $FINALFILE ]; then
      rm -f "${BAMSTOMERGE[@]}"
      rm -f $(echo "${BAMSTOMERGE[@]}" | sed 's/bam$/bai/')
    fi

    if [ -f $FINALFILE ]; then
      rm -f $DDFILE
      rm -f $(echo $DDFILE | sed 's/bam$/bai/')
    fi
  fi #FINALFILE created or already exists
  FQ=$(basename $FASTQEXT .gz)
  UMAP=$UNMAPPED/${STRAIN}.$FQ
  UMAPSINGLE=$UNMAPPED/${STRAIN}_single.$FQ
  #echo "$UMAP $UMAPSINGLE $FQ"

  if [ ! -f $UMAP.gz ]; then
    module load BBMap
    samtools fastq -f 4 --threads $CPU -N -s $UMAPSINGLE -o $UMAP $FINALFILE
    pigz $UMAPSINGLE
    repair.sh in=$UMAP out=$UMAP.gz
    unlink $UMAP
  fi
done
