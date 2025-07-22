#!/usr/bin/bash -l

echo "Strain,Pop" > PopAssigned.csv
grep immitis samples.csv | cut -d, -f2 | perl -p -e 's/^(\S+)/$1,1/' >> PopAssigned.csv
grep posadasii samples.csv | cut -d, -f2 | perl -p -e 's/^(\S+)/$1,2/' >> PopAssigned.csv

