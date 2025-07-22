#!/usr/bin/bash -l
#SBATCH -c 8 -N 1 -n 1 --mem 48gb --out logs/structure.log

module load plink
module load yq
module load faststructure
module load bcftools

OUT=structure

mkdir -p $OUT


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

for POPNAME in $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//')
do
    # create a filtered VCF containing only invariant sites
    for TYPE in SNP
    do
	for WINDOW in prune_window100 prune_ld
    	do
	    BCF=$FINALVCF/$PREFIX.$POPNAME.$TYPE.$WINDOW.bcf
	    if [ ! -f $BCF ]; then 
		echo "no pruned BCF file did you both of step 5 (05_prune_bcftools.sh) (05_pruneld_bcftools.sh)?"
		continue
	    fi
	    plink --bcf $BCF --const-fid --allow-extra-chr  --vcf-idspace-to _ --keep-allele-order --make-bed \
		  --out $OUT/$PREFIX.$POPNAME.$TYPE.$WINDOW
	    parallel -j 8 structure.py -K {} --seed 121 --input=$OUT/$PREFIX.$POPNAME.$TYPE.$WINDOW \
		     --output=$OUT/$PREFIX.$POPNAME.$TYPE.$WINDOW ::: $(seq 8)
	    chooseK.py --input=$OUT/$PREFIX.$POPNAME.$TYPE.prune_w100 > $OUT/Kchoice_$PREFIX.$POPNAME.$TYPE.$WINDOW.info
    	    pigz -kf $OUT/$PREFIX.$POPNAME.$TYPE.$WINDOW.*.meanQ
    	    Rscript scripts/plot_distruct.R $PREFIX.$POPNAME.$TYPE.$WINDOW $OUT $POPNAME
    	done
    done
done
