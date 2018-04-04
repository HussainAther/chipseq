from snakemake.shell import shell

log = snakemake.log_fmt_shell()

shell(
	"(bwa mem {snakemake.params} -t {snakemake.threads} "
	"{snakemake.input.ref} {snakemake.input.sample} "
	"| samtools view -Sbh -o {snakemake.output[0]} -) {log}")
