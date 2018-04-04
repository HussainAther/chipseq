from snakemake import shell
shell('samtools index {snakemake.input}')
shell('bamCoverage --normalizeUsingRPKM -b {snakemake.input} -o {snakemake.output}')
