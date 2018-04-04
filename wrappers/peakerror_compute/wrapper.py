from snakemake import shell

shell("""awk '{{OFS="\t"; print $1,$2,$3}}' {snakemake.input.macs2_q} | grep '^{snakemake.params.chr}.*$' > {snakemake.output.cut}""")

shell('source activate peakerror; Rscript peakerror/PeakError_compute.R {snakemake.output.cut} {snakemake.params.labels} > {snakemake.output.errors}')
