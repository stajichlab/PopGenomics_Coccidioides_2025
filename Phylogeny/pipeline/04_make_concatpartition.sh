#!/usr/bin/bash -l
#SBATCH -p short -c 2 --mem 16gb --out logs/make_mfa_cds.%A.log
CPU=${SLURM_CPUS_ON_NODE}
if [ -z $CPU ]; then
    CPU=1
fi
CPURUN=96

TIME=48:00:00

module load phykit
PREFIX=Cocci_cds

marker=ascomycota_odb10
shortmarker=$(basename $marker _odb10)
TYPE=cds
OUTDIR=results
TREEFOLDER=$OUTDIR/tree_${TYPE}_${shortmarker}

USERTREE=$TREEFOLDER/final_tree.nw
FILTERDIR=msa_filter_${TYPE}_${shortmarker}
STEM=$(ls input/*.fa | wc -l | awk '{print $1}')taxa_${shortmarker}
if [ ! -d $FILTERDIR ]; then
	echo "no $FILTERDIR folder"
	continue
fi
pushd $OUTDIR/${FILTERDIR}
ls *.mfa > filenames
mkdir -p ../${FILTERDIR}-buildtree
phykit create_concat -a filenames -p ../$FILTERDIR-buildtree/${PREFIX}.${STEM}
popd
pushd $OUTDIR/$FILTERDIR-buildtree
perl -i -p -e 's/AUTO/DNA/' ${PREFIX}.${STEM}.partition
sbatch  --time $TIME -c $CPURUN --mem 128gb -J modeltest$TYPE --out modeltest-${TYPE}.%A.log --wrap "hostname; module load modeltest-ng; modeltest-ng -i ${PREFIX}.${STEM}.fa -q ${PREFIX}.${STEM}.partition --processes $CPURUN -T raxml -d nt -t ml -T raxml"
popd
