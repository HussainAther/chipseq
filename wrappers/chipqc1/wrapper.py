from snakemake import shell
import csv

#Run chipqc on the input using the sample-wise approach
shell('module load R; Rscript ChIPQCrds.R -i={snakemake.input.reads} -p={snakemake.input.spp_peaks} -o={snakemake.output.spp}')
shell('module load R; Rscript ChIPQCrds.R -i={snakemake.input.reads} -p={snakemake.input.spp2_peaks} -o={snakemake.output.spp2}')
shell('module load R; Rscript ChIPQCrds.R -i={snakemake.input.reads} -p={snakemake.input.macs2_peaks} -o={snakemake.output.macs2}')

