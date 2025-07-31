#!/bin/bash -l
#SBATCH -N 1 -n 24 --mem 64gb --out logs/AAFTF_mito.%a.log -a 1-510

# requires AAFTF 0.3.1 or later for full support of fastp options used

MEM=64
CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

module load workspace/scratch
MITOREF=lib/MT_WA211.fa
FASTQ=input
SAMPFILE=samples.csv
ASM=asm/mito
TRIMFOLDER=$SCRATCH
WORKDIR=$SCRATCH
PHYLUM=Ascomycota
mkdir -p $ASM $TRIMFOLDER
if [ -z $CPU ]; then
    CPU=1
fi
IFS=, # set the delimiter to be ,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read RUNACC STRAIN BIOSAMPLE CENTER EXPERIMENT PROJECT ORGANISM FILEBASE NOTES LOCUSTAG
do
    module load fastp
    module load AAFTF
    ASMFILE=$ASM/${STRAIN}.mitochondria.fasta
    CLEANDUP=$ASM/${STRAIN}.rmdup.fasta
    POLISH=$ASM/${STRAIN}.polished.fasta
    SORTED=$ASM/${STRAIN}.sorted.fasta
    STATS=$ASM/${STRAIN}.sorted.stats.txt
    LEFTIN=$(ls $FASTQ/${RUNACC}${FILEBASE} | sed -n 1p)
    RIGHTIN=$(ls $FASTQ/${RUNACC}${FILEBASE} | sed -n 2p)
    BASE=$STRAIN
    
    LEFTTRIM=$WORKDIR/${BASE}_mito_1P.fastq.gz
    RIGHTTRIM=$WORKDIR/${BASE}_mito_2P.fastq.gz
    LEFT=$TRIMFOLDER/${BASE}_mito_filtered_1.fastq.gz
    RIGHT=$TRIMFOLDER/${BASE}_mito_filtered_2.fastq.gz
    echo "$RUNACC $STRAIN $LEFTIN $RIGHTIN"
    if [ ! -f $ASMFILE ]; then # can skip we already have made an assembly
	if [ ! -f $LEFT ]; then
	    if [ ! -f $LEFTTRIM ]; then
		AAFTF trim --method fastp --dedup --memory $MEM --left $LEFTIN --right $RIGHTIN -c $CPU -o $WORKDIR/${BASE}_fastp
#		AAFTF trim --method fastp --cutright -c $CPU --memory $MEM --left $WORKDIR/${BASE}_fastp_1P.fastq.gz --right $WORKDIR/${BASE}_fastp_2P.fastq.gz -o $WORKDIR/${BASE}_fastp2
		AAFTF trim --method bbduk -c $CPU --memory $MEM --left $WORKDIR/${BASE}_fastp_1P.fastq.gz --right $WORKDIR/${BASE}_fastp_2P.fastq.gz -o $WORKDIR/${BASE}_mito
	    fi
	    AAFTF filter -c $CPU --memory $MEM -o $TRIMFOLDER/${BASE}_mito --left $LEFTTRIM --right $RIGHTTRIM --aligner bbduk
	    if [ -f $LEFT ]; then
		rm -f $LEFTTRIM $RIGHTTRIM $WORKDIR/${BASE}_fastp* 
	    else
		echo "did not create left file ($LEFT $RIGHT)"
		exit
	    fi
	fi
	
	AAFTF mito --left $LEFT --right $RIGHT -o $ASMFILE -w $WORKDIR/mito_${STRAIN} --reference $MITOREF
	
	if [ ! -f $ASMFILE ]; then
	    echo "mito must have failed, exiting"
	    exit
	fi
    fi
    
    if [ ! -f $CLEANDUP ]; then
    	AAFTF rmdup -i $ASMFILE -o $CLEANDUP -c $CPU -m 500
    fi
    
    if [ ! -f $POLISH ]; then
	#AAFTF polish --method polca -i $CLEANDUP -o $POLISH -c $CPU --left $LEFT  --right $RIGHT --mem $MEM    	

	AAFTF polish --method pilon --iterations 2 -i $CLEANDUP -o $POLISH -c $CPU --left $LEFT  --right $RIGHT --mem $MEM    	
    fi
    
    if [ ! -f $POLISH ]; then
    	echo "Error running polisg, did not create file. Exiting"
    	exit
    fi
    
    if [ ! -f $SORTED ]; then
	 AAFTF sort -i $POLISH -o $SORTED
    fi
    
    if [ ! -f $STATS ]; then
	AAFTF assess -i $SORTED -r $STATS
    fi
done
