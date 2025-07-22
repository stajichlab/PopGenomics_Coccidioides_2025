#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 64 --mem 64gb --out logs/18_pairiwise_SNPs.log

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

pairwise_count() {
    strain1=$1
    strain2=$2
    TYPE=$3
    IN=$4
    if [[ "$strain1" != "$strain2" ]]
    then
        bcftools view -s $strain1,$strain2 -Ou $IN |
            bcftools +fill-tags - -Ou -- -t all | \
            bcftools view --exclude='AC=0' -f 'PASS,.' -Ou | \
            bcftools query -f '%CHROM\t%POS\t%REF[\t%TGT]\n' -o $SCRATCH/$strain1-$strain2.$TYPE.tsv
        M=$(./scripts/count_pairwise_vcftab.py --input $SCRATCH/$strain1-$strain2.$TYPE.tsv)
        echo -e "$strain1\t$strain2\t\t$M"
        rm -f $SCRATCH/$strain1-$strain2.$TYPE.tsv 
    fi
}

export -f pairwise_count

OUTDIR=reports/pairwise_strain_compare
mkdir -p $OUTDIR
for POPNAME in $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//')
do
    for TYPE in SNP INDEL
    do
        OUT=$OUTDIR/$PREFIX.$POPNAME.$TYPE.pairwise_count.tsv
	BCFILE=$PREFIX.$POPNAME.$TYPE.bcf
        IN=$FINALVCF/$BCFILE
	rsync -a $IN $SCRATCH
        echo -e "STRAIN1\tSTRAIN2\tCOUNT" > $OUT
	bcftools query -l $SCRATCH/$BCFILE > $SCRATCH/names.txt
        parallel -j $CPU pairwise_count ::: $(cat $SCRATCH/names.txt) ::: $(cat $SCRATCH/names.txt) ::: $TYPE ::: $SCRATCH/$BCFILE >> $OUT
    done
done
