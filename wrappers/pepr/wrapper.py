import tempfile
from snakemake import shell

script = """
PePr -c {snakemake.input.ip} -i {snakemake.input.inp} -f bam
""".format(**locals())

script_filename = snakemake.output[0]
fout=open(script_filename,"w")
fout.write(script)
fout.close()
shell("./{script_filename} &> {snakemake.output}.log')
