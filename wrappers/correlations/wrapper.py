from snakemake import shell

# multiBgwigSummary compares avergae scores for bigwig files which are then plotted using plotCorrelation.
shell('module load R; multiBigwigSummary bins -b {snakemake.input.ip} {snakemake.input.inp} -out {snakemake.output.mbs_spp} --BED {snakemake.input.spp_region}')
shell('module load R; multiBigwigSummary bins -b {snakemake.input.ip} {snakemake.input.inp} -out {snakemake.output.mbs_macs2} --BED {snakemake.input.macs2_region}')
shell('module load R; plotCorrelation --corData {snakemake.output.mbs_spp} -o {snakemake.output.correlations_heat_spp} -c pearson -p heatmap')
shell('module load R; plotCorrelation --corData {snakemake.output.mbs_spp} -o {snakemake.output.correlations_scat_spp} -c pearson -p scatterplot')

shell('module load R; plotCorrelation --corData {snakemake.output.mbs_macs2} -o {snakemake.output.correlations_heat_macs2} -c pearson -p heatmap')
shell('module load R; plotCorrelation --corData {snakemake.output.mbs_macs2} -o {snakemake.output.correlations_scat_macs2} -c pearson -p scatterplot')
