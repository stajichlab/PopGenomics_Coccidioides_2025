#!/usr/bin/bash -l 
#SBATCH -p short --mem 2gb --out logs/00_index.log
module load samtools
module load bwa-mem2
if [ -f config.txt ]; then
	source config.txt
fi

FASTAFILE=$REFGENOME

if [[ ! -f $FASTAFILE.fai || $FASTAFILE -nt $FASTAFILE.fai ]]; then
	samtools faidx $FASTAFILE
fi
if [[ ! -f $FASTAFILE.0123 || $FASTAFILE -nt $FASTAFILE.0123 ]]; then
	bwa-mem2 index $FASTAFILE
fi

DICT=$(dirname $FASTAFILE)/$(basename $FASTAFILE .fasta)".dict"
if [[ ! -f $DICT || $FASTAFILE -nt $DICT ]]; then
	rm -f $DICT
	samtools dict $FASTAFILE > $DICT
	ln -s $DICT $FASTAFILE.dict 
fi
#grep ">" $FASTAFILE | perl -p -e 's/>(scaffold_(\d+))/>$1,$2/' > $(dirname $FASTAFILE)/chrom_nums.csv
