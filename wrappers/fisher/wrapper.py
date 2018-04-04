from snakemake import shell

# Perform fisher's exact test on the bed files
shell('bedtools fisher -a {snakemake.input.macs2_sorted} -b {snakemake.input.spp_sorted} -g data/dm6.chr.size > {snakemake.output.fisher}')


