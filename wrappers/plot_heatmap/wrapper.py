from snakemake import shell
import os
import sys

if os.stat(snakemake.input.gzipped_matrix).st_size == 0:
    open(snakemake.output.pdf,'a').close()
else:
    shell(
    "plotHeatmap -m {snakemake.input.gzipped_matrix} " \
    "-T {snakemake.params.subprogram}_{snakemake.params.kind}_{snakemake.params.sample} " \
    "-out {snakemake.output.pdf} " \
    )
