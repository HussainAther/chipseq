from snakemake.shell import shell

shell(
    """
    bowtie2-build {snakemake.input} {snakemake.params.prefix}
    """
)

