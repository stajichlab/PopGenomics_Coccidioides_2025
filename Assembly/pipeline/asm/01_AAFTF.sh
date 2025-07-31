#!/bin/bash -l
#SBATCH -N 1 -n 1 -c 24 --mem 96gb --out logs/AAFTF.%a.log

# requires AAFTF 0.3.1 or later for full support of fastp options used
module load singularity
MEM=96
CPU=$SLURM_CPUS_ON_NODE
N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

FASTQ=input
SAMPFILE=samples.csv
ASM=asm/AAFTF
WORKDIR=$SCRATCH
#WORKDIR=working
PHYLUM=Ascomycota
mkdir -p $ASM $WORKDIR
if [ -z $CPU ]; then
    CPU=1
fi
IFS=, # set the delimiter to be ,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read RUNACC STRAIN BIOSAMPLE CENTER EXPERIMENT PROJECT ORGANISM FILEBASE NOTES LOCUSTAG
do
    ID=$STRAIN
    BASE=$STRAIN
    module load AAFTF
    module load fastp
    ASMFILE=$ASM/${ID}.spades.fasta
    VECCLEAN=$ASM/${ID}.vecscreen.fasta
    PURGE=$ASM/${ID}.sourpurge.fasta
    CLEANDUP=$ASM/${ID}.rmdup.fasta
    POLISHED=$ASM/${ID}.polished.fasta
    SORTED=$ASM/${ID}.sorted.fasta
    STATS=$ASM/${ID}.sorted.stats.txt
    LEFTIN=$(ls $FASTQ/${RUNACC}${FILEBASE} | sed -n 1p)
    RIGHTIN=$(ls $FASTQ/${RUNACC}${FILEBASE} | sed -n 2p)
    BASE=$STRAIN

    if [[ -s $SORTED ]]; then
	OLDER=1
	if [[ $LEFTIN -nt $SORTED ]]; then
	    OLDER=0
	fi
	if [[ $OLDER == 1 ]]; then
	    echo "skipping $ID $STRAIN --> $SORTED already exists and all FASTQ are older"
	    exit
	fi
    fi

    LEFTTRIM=$WORKDIR/${ID}_1P.fastq.gz
    RIGHTTRIM=$WORKDIR/${ID}_2P.fastq.gz
    MERGETRIM=$WORKDIR/${ID}_fastp_MG.fastq.gz

    # these are final processed files for assembly
    LEFT=$WORKDIR/${ID}_filtered_1.fastq.gz
    RIGHT=$WORKDIR/${ID}_filtered_2.fastq.gz
    MERGED=$WORKDIR/${ID}_filtered_U.fastq.gz

    echo "$BASE $ID $STRAIN"
    echo "$LEFTIN $RIGHTIN $LEFTTRIM $RIGHTTRIM"

    
    if [ ! -f $LEFT ]; then
		if [ ! -f $LEFTTRIM ]; then
			AAFTF trim --method fastp --dedup --merge --memory $MEM --left $LEFTIN --right $RIGHTIN -c $CPU -o $WORKDIR/${ID}_fastp
			#AAFTF trim --method fastp --cutright -c $CPU --memory $MEM --left $WORKDIR/${ID}_fastp_1P.fastq.gz --right $WORKDIR/${ID}_fastp_2P.fastq.gz -o $WORKDIR/${ID}_fastp2
			AAFTF trim --method bbduk -c $CPU --memory $MEM --left $WORKDIR/${ID}_fastp_1P.fastq.gz --right $WORKDIR/${ID}_fastp_2P.fastq.gz -o $WORKDIR/${ID}
		fi
		AAFTF filter -c $CPU --memory $MEM -o $WORKDIR/${ID} --left $LEFTTRIM --right $RIGHTTRIM --aligner bbduk
		AAFTF filter -c $CPU --memory $MEM -o $WORKDIR/${ID} --left $MERGETRIM --aligner bbduk
		if [ -f $LEFT ]; then
			rm -f $LEFTTRIM $RIGHTTRIM $WORKDIR/${ID}_fastp* 
			echo "found $LEFT"
		else
			echo "did not create left file ($LEFT $RIGHT)"
			exit
		fi
    fi
    if [ ! -f $ASMFILE ]; then # can skip we already have made an assembly
	AAFTF assemble -c $CPU --left $LEFT --right $RIGHT --merged $MERGED --memory $MEM \
	      -o $ASMFILE -w $WORKDIR/spades_${ID} --isolate --method spades -v
	
	#if [ -s $ASMFILE ]; then
	#    rm -rf $WORKDIR/spades_${ID}/K?? $WORKDIR/spades_${ID}/tmp $WORKDIR/spades_${ID}/K???
	#    rm -rf $WORKDIR/spades_${ID}
	#fi
	
	if [ ! -f $ASMFILE ]; then
	    echo "SPADES must have failed, exiting"
	    tail -n 100 $WORKDIR/spades_${ID}/spades.log
	    exit
	fi
    fi
    
    if [[ ! -f $VECCLEAN && ! -f $VECCLEAN.gz ]]; then
	AAFTF fcs_screen -i $ASMFILE -o $VECCLEAN.fcs_screen
	AAFTF vecscreen -i $VECCLEAN.fcs_screen -c $CPU -o $VECCLEAN
    fi
    if [[ ! -f $PURGE && ! -f $PURGE.gz ]]; then
	AAFTF sourpurge -i $VECCLEAN -o $PURGE -c $CPU --phylum $PHYLUM 
	# let's not remove based on coverage for now this maybe too agressive
	#--left $LEFT --right $RIGHT
	pigz $VECCLEAN
	pigz $VECCLEAN.fcs_screen 
    fi
    
    if [[ ! -f $CLEANDUP && ! -f $CLEANDUP.gz ]]; then
    	AAFTF rmdup -i $PURGE -o $CLEANDUP -c $CPU -m 500
	pigz $PURGE 
    fi
    
    if [[ ! -f $POLISHED && ! -f $POLISHED.gz ]]; then
	if [ ! -f $CLEANDUP ]; then
	    gunzip $CLEANDUP
	fi
    	AAFTF polish --method polca -i $CLEANDUP -o $POLISHED -c $CPU --left $LEFT  --right $RIGHT --mem $MEM
	pigz $CLEANDUP
    fi
    
    if [[ ! -f $POLISHED && ! -f $POLISHED.gz ]]; then
	echo "Error running polishing (polca), did not create file. Exiting"
    	exit
    fi
    
    if [ ! -f $SORTED ]; then
	AAFTF sort -i $POLISHED -o $SORTED
	pigz $POLISHED
    fi
    
    if [ ! -f $STATS ]; then
	AAFTF assess -i $SORTED -r $STATS
    fi
done
