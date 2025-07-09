#!/usr/bin/bash -l
#SBATCH --mem 64gb -N 1 -n 4 --out logs/prune_bcftools.log -p epyc

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
      OUT=$FINALVCF/$PREFIX.$POPNAME.$TYPE.prune_window100.bcf
      QC=$FINALVCF/$PREFIX.$POPNAME.$TYPE.prune_window100.stats
      if [ ! -s $OUT ];  then
	  bcftools filter -sLowQual -g3 -G10 -Ob -e 'AC==0 || AC==AN || F_MISSING > 0.02' $IN | bcftools +prune -w 100bp -n 1 -N 1st -Ob -o $OUT
     	  tabix $OUT
      fi
      if [[ ! -s $QC || $OUT -nt $QC ]]; then
	  bcftools stats $OUT > $QC
      fi 
  done
 done
