from snakemake import shell

shell('sort -k 1,1 -k2,2n {snakemake.input.beds} > {snakemake.output.sorted}')
