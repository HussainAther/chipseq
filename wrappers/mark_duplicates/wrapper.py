from snakemake.shell import shell

shell("picard MarkDuplicates {snakemake.params} INPUT={snakemake.input} "
      "OUTPUT={snakemake.output.bam} METRICS_FILE={snakemake.output.metrics} "
      "&> {snakemake.log}")
