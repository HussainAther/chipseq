from Bio import SeqIO
import gzip

chroms_to_keep=["chr2L","chr2R","chr3L","chr3R","chr4","chrX","chrY","mitochondrion_genome"]
def example_chrom_func(x):
    return x.name
def example_modifier(x):
    if x.name=="mitochondrion_genome":
        orig="M"
    else:
        orig=x.name
    new=orig
    x.id=new
    x.name=new
    x.description=x.description.replace(orig,new)
    return x
def fasta_filter(infile,outfile,chrom_func,modifier=None):
    def _iter():
        parser=SeqIO.parse(open(infile,"rt"),"fasta")
        for rec in parser:
            if chrom_func(rec) not in chroms_to_keep:
                continue
            if modifier:
                rec=modifier(rec)
            yield rec
    with open(outfile, "wt") as fout:
        SeqIO.write(_iter(),fout,"fasta")
fasta_filter('/data/Lei/Hussain/ChipSeq/data/dm6.fa', '/data/Lei/Hussain/ChipSeq/data/dm6.fixed.fa', chrom_func=example_chrom_func, modifier=example_modifier)
