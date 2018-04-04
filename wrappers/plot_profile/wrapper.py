from snakemake import shell
import os
import sys
import pybedtools

if os.stat(snakemake.input.gzipped_matrix).st_size == 0:
    open(snakemake.output.pdf, 'a').close()
    open(snakemake.output.tsv, 'a').close()
else:
    shell(
    "plotProfile -m {snakemake.input.gzipped_matrix} " \
    "-out {snakemake.output.pdf} " \
    "--outFileNameData {snakemake.output.tsv} " \
    )
