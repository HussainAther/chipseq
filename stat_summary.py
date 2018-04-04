from __future__ import division
from matplotlib.colors import LogNorm
from matplotlib.mlab import bivariate_normal
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
sample_id=[]

# This file stat_summary.py will look at the output 
# from bedtools intersect and create a tsv table of
# how many reads were overlapping, non-overlapping, etc.
#
# First we read in the input and ip files of the samples.
metadata=pd.read_csv('config/metadata.tsv',header=0,sep='\t')
metadata=metadata+"_R1"
temp_dict=pd.DataFrame.to_dict(metadata,orient='list')
d=dict(zip(temp_dict['ip'],temp_dict['input']))
d["SRR567586"]="SRR567585"

# These dictionaries will be used for storing the data 

stats = {"spp peaks" : ["spp peaks"], "macs2 peaks" : ["macs2 peaks"], "spp overlap" : ["spp overlap"], "macs2 overlap" : ["macs2 overlap"], "spp only" : ["spp only"], "macs2 only" : ["macs2 only"]}
percent_stats = {"% spp overlap" : ["% spp overlap"], "% macs2 overlap" : ["% macs2 overlap"], "% spp only" : ["% spp only"], "% macs2 only" : ["% macs2 only"]}

# This function 'file_len' returns the line count of 
# the files. Use it to determine things like the number
# of overlaps between peak-callers.
def file_len(fname): 
    with open(fname) as f:
        i=-1
        for i, l in enumerate(f):
            pass
    return int((i + 1))

# This for loop will loop through the list of sample
# names and collect the bedtools intersect information
# from them.
for i in sorted(d.keys()):
    sample_id.append(str(i))
    stats["spp peaks"].append(
        str((file_len("peak_out/spp/"
        + str(i) + "/" + str(i)
        + ".bed.narrowPeak"))))
    stats["macs2 peaks"].append(
        str((file_len("peak_out/macs2/"
        + str(i) + "/" + str(i)
        + "_peaks.narrowPeak"))))
    stats["spp overlap"].append(
        str((file_len("bedtools_out/overlap_"
        + str(i) + "_spp_macs2.txt"))))
    stats["macs2 overlap"].append(
        str((file_len("bedtools_out/overlap_"
        + str(i) + "_macs2_spp.txt"))))
    stats["spp only"].append(
        str((file_len("bedtools_out/spp_only_"
        + str(i) + "_spp_macs2.txt"))))
    stats["macs2 only"].append(
        str((file_len("bedtools_out/macs2_only_"
        + str(i) + "_spp_macs2.txt"))))

# This does the same but calculates the percentage statistics.
# These percentage overlaps are appended to the lists
#

for i in range(1,len(d.keys())+1):
    total = int(stats["spp peaks"][i]) + int(stats["macs2 peaks"][i])
    percent_stats["% spp overlap"].append(str(int(stats["spp overlap"][i]) / total))
    percent_stats["% macs2 overlap"].append(str(int(stats["macs2 overlap"][i]) / total))
    percent_stats["% spp only"].append(str(int(stats["spp only"][i]) / total))
    percent_stats["% macs2 only"].append(str(int(stats["macs2 only"][i]) / total))

# Now create the tsv file with those statistics

tsv ="\n".join(["\t".join(k) for k in
    [sample_id,stats["spp peaks"], stats["macs2 peaks"],
    stats["spp overlap"],stats["macs2 overlap"],
    stats["spp only"], stats["macs2 only"]]]) 

percent_tsv="\n".join(["\t".join(k) for k in
    [sample_id,percent_stats["% spp overlap"], percent_stats["% macs2 overlap"], 
    percent_stats["% spp only"], percent_stats["% macs2 only"]]])

with open("bedtools_out/bedtools_stats.tsv", "w") as file:
    file.write(tsv)

with open("bedtools_out/bedtools_percent_stats.tsv", "w") as file:
    file.write(percent_tsv)

# Now create the heatmap or clustermap

array=pd.DataFrame.from_csv('bedtools_out/bedtools_stats.tsv', sep='\t')
fig=plt.figure(figsize=(12,5))
overlap_heatmap=sns.heatmap(array,square=True)
fig.savefig('bedtools_out/heatmap.png')

percent_array=pd.DataFrame.from_csv('bedtools_out/bedtools_percent_stats.tsv', sep='\t')
fig=plt.figure(figsize=(12,5))
overlap_heatmap=sns.heatmap(percent_array,square=True)
fig.savefig('bedtools_out/percent_heatrmap.png')

array=pd.DataFrame.from_csv('bedtools_out/bedtools_stats.tsv', sep='\t').fillna(0)
fig=plt.figure(figsize=(12,5))
clustermap=sns.clustermap(array)
clustermap.savefig('bedtools_out/clustermap.png')

percent_array=pd.DataFrame.from_csv('bedtools_out/bedtools_percent_stats.tsv', sep='\t').fillna(0)
fig = plt.figure(figsize = (12,5))
percent_clustermap=sns.clustermap(percent_array)
percent_clustermap.savefig('bedtools_out/percent_clustermap.png')
