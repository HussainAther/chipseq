from matplotlib.colors import LogNorm
from matplotlib.mlab import bivariate_normal
from snakemake import shell
import matplotlib
matplotlib.use('Agg')
import numpy
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import scipy
import pylab
sns.set()
d_id=[]

# This file will look at the output 
# from bedtools intersect and create a tsv table of
# how many reads were overlapping, non-overlapping, etc.
#
# First we read in the input and ip files of the ds.
metadata=pd.read_csv('config/metadata.tsv',header=0,sep='\t')
metadata=metadata+"_R1"
temp_dict=pd.DataFrame.to_dict(metadata,orient='list')
d = dict(zip(temp_dict['ip'],[]*len(temp_dict['ip'])))

d = {"lm32_kc_suhw_R1" : [], "sjl_cl8_suhw_R1" : []}

# These dictionaries will be used for storing the data 
kinds = [
    'spp',
    'macs2',
    'spp2',
    'macs2_only',
    'spp_only',
    'spp2_only',
    'macs2_found-in-all',
    'spp_found-in-all',
    'spp2_found-in-all',
    'macs2-not-spp',
    'macs2-not-spp2',
    'spp-not-macs2',
    'spp-not-spp2',
    'spp2-not-macs2',
    'spp2-not-spp'
]

stats = {
"spp" : d, 
"spp2": d, 
"macs2" : d, 
"spp-not-spp2" : d, 
"spp2-not-spp": d, 
"macs2-not-spp" : d,
"macs2-not-spp2" : d,
"spp-not-macs2" : d,
"spp2-not-macs2" : d,
"spp_only" : d, 
"spp2_only": d, 
"macs2_only" : d,
"macs2_found-in-all" : d,
"spp_found-in-all" : d,
"spp2_found-in-all" : d
}

stats2 = stats

f = open(snakemake.input.profile_ref)
f2 = open(snakemake.input.profile_sca)

# This iterates through the files and collects the coverage values
# for the genes at each position from the TSS

for i, line in enumerate(f):
    if i == 3:
        for i in line.split("\t")[2:]:
            stats[{snakemake.wildcards.kind}][{snakemake.wildcards.sample}].append(i)

for i, line in enumerate(f2):
    if i == 3:
        for i in line.split("\t")[2:]:
            stats2[{snakemake.wildcards.kind}][{snakemake.wildcards.sample}].append(i)

output = ""

for kind in stats:
    for sample in kind:
       output += stats[kind][sample] + "\t"
    output += "\n"

f=open(snakemake.output.ref,"w")
f.write(str(output))
f.close()

output = ""

for kind in stats:
    for sample in kind:
       output += stats2[kind][sample] + "\t"
    output += "\n"

f=open(snakemake.output.sca,"w")
f.write(str(output))
f.close()
