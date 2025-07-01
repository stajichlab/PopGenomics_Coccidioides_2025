#!/usr/bin/bash -l

echo -e "Populations:\n\n  All:" > population_sets.yaml
IFS=,
tail -n +2 samples.csv | while read RUN STRAIN OTHER FILES
do
	echo "    - $STRAIN"
done >> population_sets.yaml
