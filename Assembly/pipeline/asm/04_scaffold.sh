#!/usr/bin/bash -l
#SBATCH -p short -c 64 --mem 64gb --out logs/ragtag.log

module load ragtag
module load parallel
CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
    CPU=1
fi
CIREF=../Genotyping/genome/FungiDB-68_CimmitisRS_Genome.fasta
CPREF=../Genotyping/genome/FungiDB-68_CposadasiiSilveira2022_Genome.fasta
OUTDIR=asm/scaffold
INDIR=genomes

mkdir -p "$OUTDIR"

run_ragtag() {
    ref=$1
    target=$2
    outdir=$3

    if [ ! -f "$ref" ]; then
        echo "Ref genome $ref does not exist."
        return 1
    fi

    if [ ! -f "$target" ]; then
        echo "Target genome file $target does not exist."
        return 1
    fi
    if [ -z "$outdir" ]; then
        echo "Output file not specified."
        return 1
    fi
    if [ -d "$outdir" ]; then
        echo "Output folder $outdir already exists. Skipping."
        return 0
    fi
    ragtag.py scaffold -u -t 2 "$ref" "$target" -o "$outdir"
}

export -f run_ragtag
parallel -j 4 run_ragtag $CIREF $INDIR/{}.AAFTF.fasta "$OUTDIR"/{} ::: $(tail -n +2 samples.csv | grep "Coccidioides immitis" | cut -d, -f2)
parallel -j 4 run_ragtag $CPREF $INDIR/{}.AAFTF.fasta "$OUTDIR"/{} ::: $(tail -n +2 samples.csv | grep "Coccidioides posadasii" | cut -d, -f2)