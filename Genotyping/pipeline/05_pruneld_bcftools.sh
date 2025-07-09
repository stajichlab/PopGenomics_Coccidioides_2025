#!/usr/bin/bash -l
#SBATCH --mem 64gb -N 1 -n 4 --out logs/pruneLD_bcftools.log -p epyc

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
if [[ -z $POPYAML || ! -s $POPYAML ]]; then
    echo "Cannot find \$POPYAML variable - set in config.txt"
    exit
fi

for POPNAME in $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//')
do
  for TYPE in SNP INDEL
  do
      IN=$FINALVCF/$PREFIX.$POPNAME.$TYPE.bcf
      OUT=$FINALVCF/$PREFIX.$POPNAME.$TYPE.prune_ld.bcf
      QC=$FINALVCF/$PREFIX.$POPNAME.$TYPE.prune_ld.stats
      if [ ! -s $OUT ];  then
	 bcftools +prune -m 0.6 -w 1000 -n 1 -N 1st -Ob -o $OUT $IN
     	 tabix $OUT
      fi
      if [[ ! -s $QC || $OUT -nt $QC ]]; then
	  bcftools stats $OUT > $QC
      fi 
  done
 done
