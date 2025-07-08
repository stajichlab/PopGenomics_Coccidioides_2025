#!/usr/bin/bash -l
#SBATCH --mem 64gb -N 1 -n 4 --out logs/concat_bcf.log -p short

module load bcftools
module load yq

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi
if [ -f config.txt ]; then
	source config.txt
else
	echo "need a config.txt"
fi

if [ -z $FINALVCF ]; then
	echo "Need to define FINALVCF in config.txt"
	exit
fi
if [ -z $SLICEVCF ]; then
	echo "Need to define SLICEVCF in config.txt"
	exit
fi
if [[ -z $POPYAML || ! -s $POPYAML ]]; then
	echo "Cannot find \$POPYAML variable - set in config.txt"
	exit
fi
IN=$SLICEVCF
mkdir -p $FINALVCF

for POPNAME in $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//')
do
  for TYPE in SNP INDEL
  do
     OUT=$FINALVCF/$PREFIX.$POPNAME.$TYPE.bcf
     if [ ! -s $OUT ];  then
     	bcftools concat -Oz -o $OUT --threads $CPU $IN/$POPNAME/${PREFIX}.*.bcftools.${TYPE}.bcf
     	tabix $OUT
     fi

   done
 done
