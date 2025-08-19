#!/bin/bash -l
#SBATCH -p batch -N 1 -n 4 --mem 32gb --time 5-0:0:0 --out logs/mysql.log -J mysqld
# Define program name
PROGNAME=$(basename $0)

# Load software
module load singularity
module load workspace/scratch

# Define stop mysqldb
stop_mysqldb() { singularity instance stop mysqldb; }

# Set trap to ensure mysqldb is stopped
trap "stop_mysqldb; exit 130" SIGHUP SIGINT SIGTERM

# Define error handler
error_exit()
{
    stop_mysqldb
	echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2
	exit 1
}

# Set some vars
mkdir -p $SCRATCH/db $SCRATCH/conf

export SINGULARITY_BINDPATH=$SCRATCH
export SINGULARITYENV_PASACONF=$SCRATCH/pasa-local.config.txt
cp ~/pasa-local.config.txt $SINGULARITYENV_PASACONF || error_exit "Failed to copy pasa config file"

SIF=~/bigdata/mysql/mariadb.sif

cd $SCRATCH

# Update PASA DB config
HOSTNAME=$(hostname -s)
echo "$HOSTNAME:$PORT"
sed -i "s/^MYSQLSERVER.*$/MYSQLSERVER=${HOSTNAME}:${PORT}/" ${SINGULARITYENV_PASACONF}
singularity exec --writable-tmpfs -B db/:/var/lib/mysql $SIF /usr/bin/mysqld_safe
stop_mysqldb
exit 0
