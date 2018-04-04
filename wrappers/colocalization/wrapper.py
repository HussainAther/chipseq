from matplotlib.colors import LogNorm
from matplotlib.mlab import bivariate_normal
from snakemake import shell
import matplotlib
matplotlib.use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
sns.set()
sample_id=[]

# This file will look at the output 
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

stats = {"spp seg" : ["spp seg"], "macs2 seg" : ["macs2 seg"], "spp query": ["spp query"], "macs2 query": ["macs2 query"]}

# This function 'file_len' returns the line count of 
# the files. Use it to determine things like the number
# of overlaps between peak-callers.
def file_len(fname): 
    with open(fname) as f:
        i=-1
        for i, l in enumerate(f):
            pass
    return int((i + 1))

def comma_count(fname):
    with open(fname) as f:
        commas = 0
        for i, l in enumerate(f):
            commas += l.count(",")
    return int(commas)

# This for loop will loop through the list of sample
# names and collect the bedtools intersect information
# from them

for i in sorted(d.keys()):
    sample_id.append(str(i))
    stats["spp seg"].append(
        str((comma_count("gat/"
        + str(i)
        + ".spp.nucleotide-overlap.counts"))))
    stats["macs2 seg"].append(
        str((comma_count("gat/"
        + str(i)
        + ".macs2.nucleotide-overlap.counts"))))
    stats["spp query"].append(
        str((file_len("intervalstats/"
        + str(i)
        + "_query_spp.full"))))
    stats["macs2 query"].append(
        str((file_len("intervalstats/"
        + str(i)
        + "_query_macs2.full"))))

# Now create the tsv file with those statistics

tsv ="\n".join(["\t".join(k) for k in
    [sample_id,stats["spp seg"], stats["macs2 seg"], stats["spp query"],
    stats["macs2 query"]]]) 

with open("colocalization/colocalization_stats.tsv", "w") as file:
    file.write(tsv)

tsv ="\n".join(["\t".join(k) for k in
    [sample_id,stats["spp seg"], stats["macs2 seg"]]]) 

with open("colocalization/gat_stats.tsv", "w") as file:
    file.write(tsv)

tsv ="\n".join(["\t".join(k) for k in
    [stats["spp query"],
    stats["macs2 query"]]]) 

with open("colocalization/intervalstats_stats.tsv", "w") as file:
    file.write(tsv)

# Now create the heatmap or clustermap

array=pd.DataFrame.from_csv('colocalization/colocalization_stats.tsv', sep='\t')
fig=plt.figure(figsize=(12,5))
overlap_heatmap=sns.heatmap(array,square=True)
fig.savefig(snakemake.output.colocalization)

array=pd.DataFrame.from_csv('colocalization/gat_stats.tsv', sep='\t')
fig=plt.figure(figsize=(12,5))
overlap_heatmap=sns.heatmap(array,square=True)
fig.savefig(snakemake.output.gat)

array=pd.DataFrame.from_csv('colocalization/intervalstats_stats.tsv', sep='\t')
fig=plt.figure(figsize=(12,5))
overlap_heatmap=sns.heatmap(array,square=True)
fig.savefig(snakemake.output.intervalstats)
