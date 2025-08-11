#!/usr/bin/bash -l
#SBATCH -p epyc -c 8 --mem 4gb --out logs/buildtree_raxml.%A.log

module load raxml-ng
PREFIX=Cocci_cds

marker=ascomycota_odb10
shortmarker=$(basename $marker _odb10)
TYPE=cds
OUTDIR=results
USERTREE=$(realpath $TREEFOLDER/final_tree.nw)
FILTERDIR=msa_filter_${TYPE}_${shortmarker}
STEM=$(ls input/*.fa | wc -l | awk '{print $1}')taxa_${shortmarker}

pushd $OUTDIR/$FILTERDIR-buildtree
IN=${PREFIX}.${STEM}.fa.raxml.rba
if [ ! -f $IN ]; then
    raxml-ng --parse --model ${PREFIX}.${STEM}.fa.part.aic --msa ${PREFIX}.${STEM}.fa
fi
raxml-ng --all --msa $IN --tree pars{10} --bs-trees 200 --threads auto{7} --workers auto{4}
