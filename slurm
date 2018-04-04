#!/bin/bash
#SBATCH --job-name="snakemake"
#SBATCH --partition="norm"
#SBATCH --time=24:00:00

#interpret arguments
SNAKEFILE=$1
if [ -z $2 ]; then
    EXTRA=$2
else
     EXTRA=''
fi

# Run snakemake
(source activate snakemake; \
    time snakemake \
    -p \
    -s $SNAKEFILE \
    --directory $PWD \
    -T \
    -j 999 \
    --cluster-config config/clusterconfig.yml \
    --verbose \
    --use-conda \
    --cluster 'sbatch {cluster.prefix} --cpus-per-task={threads}' \
    --latency-wait=60 \
    $EXTRA \
) > "${SNAKEFILE}.log" 2>&1

SNAKE_PID=$!

finish(){
    echo 'Stopping running snakemake job.'
    kill -SIGINT $SNAKE_PID
    exit 0
}
trap finish SIGTERM

wait $SNAKE_PID
