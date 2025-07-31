#!/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 16G -p short --out logs/stats__busco.%a.log -J busco


module load busco
export AUGUSTUS_CONFIG_PATH=$(realpath lib/augustus/3.5/config)

module load workspace/scratch

CPU=2
if [ ! -z ${SLURM_CPUS_ON_NODE} ]; then
	CPU=${SLURM_CPUS_ON_NODE}
fi

N=${SLURM_ARRAY_TASK_ID}
if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi
GENOMEFOLDER=genomes
EXT=fasta
LINEAGE=ascomycota_odb10
OUTFOLDER=BUSCO
SAMPFILE=samples.csv
SEED_SPECIES=coccidioides_immitis
mkdir -p $OUTFOLDER

IFS=, # set the delimiter to be ,
tail -n +2 $SAMPFILE | sed -n ${N}p | while read RUNACC STRAIN BIOSAMPLE CENTER EXPERIMENT PROJECT ORGANISM FILEBASE NOTES LOCUSTAG
do
    ID=$STRAIN
    for type in AAFTF
    do
		GENOMEFILE=$GENOMEFOLDER/$ID.$type.$EXT
		if [ -f $GENOMEFILE ]; then
			echo "GENOMEFILE is $GENOMEFILE"
			GENOMEFILE=$(realpath $GENOMEFILE)
			if [ -d "$OUTFOLDER/${ID}.${type}" ];  then
				echo "Already have run $ID in folder busco - do you need to delete it to rerun?"
				exit
			else
				busco -f -m genome -l $LINEAGE -c $CPU -o ${ID}.${type} --out_path ${OUTFOLDER} --offline \
				--augustus_species $SEED_SPECIES --in $GENOMEFILE --download_path $BUSCO_LINEAGES
			fi
		fi
    done
done
