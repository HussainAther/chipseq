import os
import sys
# this file determines which input and ip correspond with one another from the metadata. the metada file is the list of input and ip samples.
# usage: python separate_ip_input.py -i METADATA -o TABLE
metadata=open("%s" % str(sys.argv[2]), "r")
table=open("%s" % str(sys.argv[4]),"w")
samples={}
for i in metadata.readlines():
    if i.split()[0]!="ip":
        samples[str(i.split()[0])]=str(i.split()[1])
metadata.close()
for i in samples:
    table.write(str(i)+str("_R1")+str("    ")+str(samples[i])+str("_R1"))
    table.write("\n")
table.close()
