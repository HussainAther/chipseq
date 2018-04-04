from snakemake.shell import shell

broad_cutoff = float(snakemake.params.qval)*1.122

shell(
    'macs2 '
    'callpeak '
    '-c {snakemake.input.inp} '
    '-t {snakemake.input.ip} '
    '--bdg --SPMR '
    '-q {snakemake.params.qval} '
    '--broad --broad-cutoff {broad_cutoff} '
    '-n {snakemake.wildcards.sample} '
    '--outdir {snakemake.params.qval_prefix}'
)
shell('Rscript {snakemake.params.sample_prefix}_model.r')
#shell("ln -sf {snakemake.output.narrowPeak} {snakemake.output.bed}")
#shell("cd {snakemake.params.qval_prefix} && ln -sf $(basename {snakemake.output.narrowPeak}) $(basename {snakemake.output.bed})")

