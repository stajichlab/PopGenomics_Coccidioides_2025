#!/usr/bin/bash -l
#SBATCH -p short -N 1 -n 1 -c 96 --mem 64gb --out logs/17_strainwise_SNPs.log

# this script generates number of SNPs each strain has with the ref

module load bcftools
module load workspace/scratch
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

count_snps() {
    strain=$1
    IN=$2
    INCPU=2
    bcftools view --threads $INCPU -s $strain -Ou $IN |
        bcftools +fill-tags --threads $INCPU - -Ob -o $SCRATCH/$strain.tags.bcf  -- -t all
    # %CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT]
    count=$(bcftools query -f '%CHROM\t%POS\n' -e "AF!=1 | FILTER!='PASS'" $SCRATCH/$strain.tags.bcf | wc -l)
    echo -e "$strain\t$count"
}

export -f count_snps
OUTDIR=reports/individual_strain_compare
mkdir -p $OUTDIR
for POPNAME in $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//')
do
    for TYPE in SNP INDEL
    do
        OUT=$OUTDIR/$PREFIX.$POPNAME.$TYPE.count.tsv
        IN=$FINALVCF/$PREFIX.$POPNAME.$TYPE.bcf
        echo -e "STRAIN\tCOUNT" > $OUT
        parallel -j $CPU count_snps ::: $(bcftools query -l $IN) ::: $IN | sort >> $OUT
        rm -f $SCRATCH/*.bcf
    done
done
