from snakemake.shell import shell

shell(
    'macs2 '
    'callpeak '
    '-c {snakemake.input.inp} '
    '-t {snakemake.input.ip} '
    '-n {snakemake.wildcards.sample} '
    '--outdir {snakemake.params.sample_prefix}'
)
#shell('Rscript {snakemake.params.sample_prefix}_model.r')
#shell("ln -sf {snakemake.output.narrowPeak} {snakemake.output.bed}")
#shell("cd {snakemake.params.qval_prefix} && ln -sf $(basename {snakemake.output.narrowPeak}) $(basename {snakemake.output.bed})")
