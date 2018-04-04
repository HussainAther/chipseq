from snakemake import shell
shell('bedtools jaccard -a {snakemake.input.macs2_sorted} -b {snakemake.input.spp_sorted} > {snakemake.output.jaccard}')
