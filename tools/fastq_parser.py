from Bio import SeqIO
count=0
for rec in SeqIO.parser("data/SRR567585.fastq", "fastq"):
    count+=1
print("%i reads" % count)

good_reads=(rec for rec in \ 
    SeqIO.parse("data/SRR567585.fastq","fastq") \
    if min(rec.letter_annotations["phred_quality"])>=20)
count=seqIO.write(good_reads, "good_quality.fastq", "fastq")
print("Saved %i reads" % count)

