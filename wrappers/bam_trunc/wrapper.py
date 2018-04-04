from snakemake import shell
from pybedtools import BedTool

# truncate the output files to chromosome so no ends are out of bounds.
BedTool(snakemake.input).truncate_to_chrom('dm6').saveas(snakemake.output)


