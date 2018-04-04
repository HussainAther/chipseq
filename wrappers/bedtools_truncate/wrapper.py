from snakemake import shell
from pybedtools import BedTool

# truncate the sorted output files to chromosome so no ends are out of bounds.
spp_trunc = snakemake.output.spp_sorted_tr
spp_sorted = snakemake.input.spp_sorted
macs2_trunc = snakemake.output.macs2_sorted_tr
macs2_sorted = snakemake.input.macs2_sorted
spp2_trunc = snakemake.output.spp2_sorted_tr
spp2_sorted = snakemake.input.spp2_sorted

genome = snakemake.input.genome

BedTool(spp_sorted).truncate_to_chrom('dm6').saveas(spp_trunc)
BedTool(spp2_sorted).truncate_to_chrom('dm6').saveas(spp2_trunc)
BedTool(macs2_sorted).truncate_to_chrom('dm6').saveas(macs2_trunc)


