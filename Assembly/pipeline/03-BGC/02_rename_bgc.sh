#!/usr/bin/bash -l
#SBATCH -c 8 --mem 16gb -p short

module load biopython
# will assume there are folder named 'antismash_local' in each strain annotation folder
python scripts/gather_rename_AS_regionfiles.py annotation results/antismash_regions
