from snakemake.shell import shell

shell(
    """
    bowtie2 \
    -x {snakemake.params.prefix} \
    -U {snakemake.input.fastq}  \
    -S {snakemake.output}.sam \
    2> {snakemake.log} && \
    samtools view -Sbh {snakemake.output}.sam > {snakemake.output} && \
    rm {snakemake.output}.sam
    """)

