from snakemake import shell
import os
import sys
import pybedtools

if os.stat(snakemake.input.regions).st_size == 0:
    open(snakemake.output.tsv_out, 'a').close()
    open(snakemake.output.gzipped_matrix, 'a').close()
else:
    shell(
    "computeMatrix {snakemake.params.subprogram} -R {snakemake.input.regions} " \
    "-S {snakemake.input.ip} {snakemake.input.inp} " \
    "--outFileNameMatrix {snakemake.output.tsv_out} " \
    "-out {snakemake.output.gzipped_matrix}"
    )
