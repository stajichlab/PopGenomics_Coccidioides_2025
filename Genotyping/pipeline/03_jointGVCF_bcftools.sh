#!/usr/bin/bash -l
#SBATCH --mem 24G --nodes 1 --ntasks 8 -p epyc -J slice.GVCFGeno --out logs/GVCFGenoGATK4_bcftools.slice_%a.%A.log  -a 1-7
hostname
MEM=24g
module unload R
module unload java
module load picard
module load gatk/4.6.1.0
module load bcftools
module load parallel
module load yq
module load workspace/scratch

source config.txt

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi
if [ -f config.txt ]; then
	source config.txt
fi
TEMPDIR=$SCRATCH
if [ ! -f $REFGENOME ]; then
    module load samtools
    samtools faidx $REFGENOME
fi
NSTART=$(perl -e "printf('%d',1 + $GVCF_INTERVAL * ($N - 1))")
NEND=$(perl -e "printf('%d',$GVCF_INTERVAL * $N)")
MAX=$(wc -l $REFGENOME.fai | awk '{print $1}')
if [ "$NSTART" -gt "$MAX" ]; then
	echo "NSTART ($NSTART) > $MAX"
	exit
fi
if [ "$NEND" -gt "$MAX" ]; then
	NEND=$MAX
fi
echo "$NSTART -> $NEND"

CPU=$SLURM_CPUS_ON_NODE
if [ ! $CPU ]; then
    CPU=2
fi
if [[ $(ls $GVCFFOLDER | grep -c -P "\.g.vcf$") -gt "0" ]]; then
   parallel -j $CPU bgzip {} ::: $GVCFFOLDER/*.g.vcf
  parallel -j $CPU tabix -f {} ::: $GVCFFOLDER/*.g.vcf.gz
fi

if [[ -z $POPYAML || ! -s $POPYAML ]]; then
	echo "Cannot find \$POPYAML variable - set in config.txt"
	exit
fi
if [ -z $SLICEVCF ]; then
	SLICEVCF=vcf_slice
fi
mkdir -p $SLICEVCF
for POPNAME in $(yq eval '.Populations | keys' $POPYAML | perl -p -e 's/^\s*\-\s*//')
do
	FILES=$(yq eval '.Populations.'$POPNAME'[]' $POPYAML | perl -p -e "s/(\S+)/-V $GVCFFOLDER\/\$1.g.vcf.gz/g"  )
	INTERVALS=$(cut -f1 $REFGENOME.fai  | sed -n "${NSTART},${NEND}p" | perl -p -e 's/(\S+)\n/--intervals $1 /g')

	mkdir -p $SLICEVCF/$POPNAME
	STEM=$SLICEVCF/$POPNAME/$PREFIX.$N
	GENOVCFOUT=$STEM.all.vcf
	SELECTSNP=$STEM.bcftools.SNP.bcf
	SELECTINDEL=$STEM.bcftools.INDEL.bcf
	echo "$STEM is stem; GENOVCFOUT=$STEM.all.vcf POPNAME=$POPNAME slice=$SLICEVCF"
	mkdir -p $TEMPDIR
	if [ ! -f $GENOVCFOUT.gz ]; then
	    if [ ! -f $GENOVCFOUT ]; then
		DB=$TEMPDIR/${GVCFFOLDER}_slice_$N
		rm -rf $DB
		gatk  --java-options "-Xmx$MEM -Xms$MEM" GenomicsDBImport --consolidate --merge-input-intervals --genomicsdb-workspace-path $DB $FILES $INTERVALS --tmp-dir $TEMPDIR --reader-threads $CPU 
		#gatk  --java-options "-Xmx$MEM -Xms$MEM" GenomicsDBImport --genomicsdb-workspace-path $DB $FILES $INTERVALS  --reader-threads $CPU
		time gatk GenotypeGVCFs --reference $REFGENOME --output $GENOVCFOUT -V gendb://$DB --tmp-dir $TEMPDIR -G StandardAnnotation -G AS_StandardAnnotation 


		ls -l $TEMPDIR
		rm -rf $DB
	    fi
	    if [ -f $GENOVCFOUT ]; then
	    	bgzip $GENOVCFOUT
	    	tabix $GENOVCFOUT.gz
	    fi
	fi
	TYPE=SNP
	echo "VCF = $STEM.$TYPE.vcf.gz"
	if [[ ! -f $STEM.$TYPE.vcf.gz ]]; then
	    gatk SelectVariants \
		-R $REFGENOME \
		--variant $GENOVCFOUT.gz \
		-O $STEM.$TYPE.vcf \
		--restrict-alleles-to BIALLELIC \
		--select-type-to-include $TYPE --create-output-variant-index false

	    bgzip $STEM.$TYPE.vcf
	    tabix $STEM.$TYPE.vcf.gz
	fi

	if [[ ! -f $SELECTSNP || $STEM.$TYPE.vcf.gz -nt $SELECTSNP ]]; then
	    bcftools filter -Ob -o $SELECTSNP -g3 -G10 -e 'QD < 2.0 || MQ < 40.0 || SOR > 3.0 || FS > 60.0 || MQRankSum < -12.5 || ReadPosRankSum < -12.5 || QUAL <10 || (AC<2 && QUAL<15)' $STEM.$TYPE.vcf.gz 
	  tabix $SELECTSNP
	  bcftools stats $SELECTSNP > $SELECTSNP.stats
	fi

	TYPE=INDEL
	if [ ! -f $STEM.$TYPE.vcf.gz ]; then
	    gatk SelectVariants \
	        -R $REFGENOME \
	        --variant $GENOVCFOUT.gz \
	        -O $STEM.$TYPE.vcf  --select-type-to-include MIXED --select-type-to-include MNP \
	        --select-type-to-include $TYPE --create-output-variant-index false
	    bgzip $STEM.$TYPE.vcf
	    tabix $STEM.$TYPE.vcf.gz
	fi


	if [[ ! -f $SELECTINDEL || $STEM.$TYPE.vcf.gz -nt $SELECTINDEL  ]]; then
	    bcftools filter -Ob -o $SELECTINDEL -sLowQual -g3 -G10 -e 'QD < 2.0 || FS > 200.0' $STEM.$TYPE.vcf.gz
	    tabix $SELECTINDEL
	    bcftools stats $SELECTINDEL > $SELECTINDEL.stats
	fi
done
