from snakemake import shell
from pybedtools import BedTool

#find the intersecting regions between the peak-callers
shell('bedtools intersect -a {snakemake.input.spp} -b {snakemake.input.macs2} -u > {snakemake.output.spp_macs2}')
shell('bedtools intersect -a {snakemake.input.spp2} -b {snakemake.input.macs2} -u > {snakemake.output.spp2_macs2}')
shell('bedtools intersect -a {snakemake.input.macs2} -b {snakemake.input.spp} -u > {snakemake.output.macs2_spp}')
shell('bedtools intersect -a {snakemake.input.macs2} -b {snakemake.input.spp2} -u > {snakemake.output.macs2_spp2}')
#shell('bedtools intersect -a {snakemake.input.macs2} -b {snakemake.input.spp} -u > {snakemake.output.macs2_pepr}')
#shell('bedtools intersect -a {snakemake.input.macs2} -b {snakemake.input.spp} -u > {snakemake.output.pepr_spp}')

#sort the output files
shell('sortBed -i {snakemake.input.spp} > {snakemake.output.spp_sorted}')
shell('sortBed -i {snakemake.input.spp2} > {snakemake.output.spp2_sorted}')
shell('sortBed -i {snakemake.input.macs2} > {snakemake.output.macs2_sorted}')
#shell('sortBed -i {snakemake.input.pepr}')

#find the non-intersecting regions between the peak-callers
shell('bedtools intersect -a {snakemake.input.spp} -b {snakemake.input.macs2} -v > {snakemake.output.spp_only}')
shell('bedtools intersect -a {snakemake.input.spp2} -b {snakemake.input.macs2} -v > {snakemake.output.spp2_only}')
shell('bedtools intersect -a {snakemake.input.macs2} -b {snakemake.input.spp} -v > {snakemake.output.macs2_only}')

