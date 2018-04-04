# ChipSeq Snakemake README
## What is this?
This is a snakemake pipeline for Chip-Seq analysis. It runs Chip-Seq commands on input sequence data and analyzes them in such a way that the peak-callers are compared. It's used with the slurm workflow management as `sbatch slurm Snakefile $METADATA` with wrappers for rule-calling. It supports parallelized jobs that can be run simultaneously. The output is analyzed using bedtools and deeptools, and the colocalized areas of the genomes are shown and compared as well. The different peak-callers (spp, macs2, etc.) are compared in greater detail. 
