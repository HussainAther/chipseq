from matplotlib.colors import LogNorm
from matplotlib.mlab import bivariate_normal
from snakemake import shell
import matplotlib
matplotlib.use('Agg')
import numpy
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
d["dc92_bg3_ctcf-56"] = "dc91_bg3_input"
d["dc94_bg3_ctcf-59"] = "dc93_bg3_input"


#sample = ["lm32_kc_suhw_R1", "sjl_cl8_suhw_R1"]

sample = []

for i in d:
    sample.append(i)

# These dictionaries will be used for storing the data 

stats = {
"spp peaks" : ["spp peaks"], 
"spp2 peaks": ["spp2 peaks"], 
"macs2 peaks" : ["macs2 peaks"], 
"spp-not-spp2" : ["spp-not-spp2"], 
"spp2-not-spp": ["spp2-not-spp"], 
"macs2-not-spp" : ["macs2-not-spp"],
"macs2-not-spp2" : ["macs2-not-spp2"],
"spp-not-macs2" : ["spp-not-macs2"],
"spp2-not-macs2" : ["spp2-not-macs2"],
"spp only" : ["spp only"], 
"spp2 only":["spp2 only"], 
"macs2 only" : ["macs2 only"],
"macs2 found in all" : ["macs2 found in all"],
"spp found in all" : ["spp found in all"],
"spp2 found in all" : ["spp2 found in all"]
}

percent_stats = {
"% spp-not-macs2" : ["% spp-not-macs2"], 
"% spp2-not-spp" : ["% spp2-not-spp"], 
"% spp-not-spp2" : ["% spp-not-spp2"], 
"% spp2-not-macs2" : ["% spp2-not-macs2"], 
"% macs2-not-spp" : ["% macs2-not-spp"], 
"% macs2-not-spp2" : ["% macs2-not-spp2"], 
"% spp_only" : ["% spp only"], 
"% spp2_only": ["% spp2 only"], 
"% macs2_only" : ["% macs2 only"],
"% macs2 found in all" : ["% macs2 found in all"],
"% spp found in all" : ["% spp found in all"],
"% spp2 found in all" : ["% spp2 found in all"]
}

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
# from them
for i in sorted(sample):
    sample_id.append(str(i))
    stats["spp peaks"].append(
        str((file_len("peak_out/spp/"
        + str(i) + "/" + str(i)
        + ".bed.narrowPeak"))))
    stats["spp2 peaks"].append(
        str((file_len("peak_out/spp2/"
        + str(i) + "/" + str(i)
        + ".bed.narrowPeak"))))
    stats["macs2 peaks"].append(
        str((file_len("peak_out/macs2/"
        + str(i) + "/" + str(i)
        + "_peaks.narrowPeak"))))
    stats["spp-not-macs2"].append(
        str((file_len("bedtools_out/spp-not-macs2/spp-not-macs2_"
        + str(i) + ".txt"))))
    stats["spp2-not-macs2"].append(
        str((file_len("bedtools_out/spp2-not-macs2/spp2-not-macs2_"
        + str(i) + ".txt"))))
    stats["spp-not-spp2"].append(
        str((file_len("bedtools_out/spp-not-spp2/spp-not-spp2_"
        + str(i) + ".txt"))))
    stats["spp2-not-spp"].append(
        str((file_len("bedtools_out/spp2-not-spp/spp2-not-spp_"
        + str(i) + ".txt"))))
    stats["macs2-not-spp"].append(
        str((file_len("bedtools_out/macs2-not-spp/macs2-not-spp_"
        + str(i) + ".txt"))))
    stats["macs2-not-spp2"].append(
        str((file_len("bedtools_out/macs2-not-spp2/macs2-not-spp2_"
        + str(i) + ".txt"))))
    stats["spp only"].append(
        str((file_len("bedtools_out/spp_only/spp_only_"
        + str(i) + ".txt"))))
    stats["macs2 only"].append(
        str((file_len("bedtools_out/macs2_only/macs2_only_"
        + str(i) + ".txt"))))
    stats["spp2 only"].append(
        str((file_len("bedtools_out/spp2_only/spp2_only_"
        + str(i) + ".txt"))))
    stats["spp found in all"].append(
        str((file_len("bedtools_out/spp_found-in-all/spp_found-in-all_"
        + str(i) + ".txt"))))
    stats["macs2 found in all"].append(
        str((file_len("bedtools_out/macs2_found-in-all/macs2_found-in-all_"
        + str(i) + ".txt"))))
    stats["spp2 found in all"].append(
        str((file_len("bedtools_out/spp2_found-in-all/spp2_found-in-all_"
        + str(i) + ".txt"))))

# This does the same but calculates the percentage statistics.
# These percentage overlaps are appended to the lists
#

for i in range(1,len(sample)+1):
    total = int(stats["spp peaks"][i]) + int(stats["macs2 peaks"][i]) + int(stats["spp2 peaks"][i])
    percent_stats["% spp-not-spp2"].append(str(int(stats["spp-not-spp2"][i]) / total))
    percent_stats["% spp2-not-spp"].append(str(int(stats["spp2-not-spp"][i]) / total)) 
    percent_stats["% spp-not-macs2"].append(str(int(stats["spp-not-macs2"][i]) / total))
    percent_stats["% spp2-not-macs2"].append(str(int(stats["spp2-not-macs2"][i]) / total)) 
    percent_stats["% macs2-not-spp"].append(str(int(stats["macs2-not-spp"][i]) / total))
    percent_stats["% macs2-not-spp2"].append(str(int(stats["macs2-not-spp2"][i]) / total))
    percent_stats["% spp_only"].append(str(int(stats["spp only"][i]) / total))
    percent_stats["% spp2_only"].append(str(int(stats["spp2 only"][i]) / total))
    percent_stats["% macs2_only"].append(str(int(stats["macs2 only"][i]) / total))
    percent_stats["% spp found in all"].append(str(int(stats["spp found in all"][i]) / total))
    percent_stats["% spp2 found in all"].append(str(int(stats["spp2 found in all"][i]) / total))
    percent_stats["% macs2 found in all"].append(str(int(stats["macs2 found in all"][i]) / total))

# Now create the tsv file with those statistics

tsv = "\n".join(["\t".join(k) for k in
    [sample_id, stats["spp peaks"], stats["macs2 peaks"], 
    stats["spp2 peaks"], stats["spp-not-macs2"], stats["spp2-not-macs2"],
    stats["macs2-not-spp"], stats["macs2-not-spp2"], stats["spp2-not-spp"],
    stats["spp-not-spp2"], stats["spp only"], 
    stats["macs2 only"], stats["spp2 only"],
    stats["spp found in all"], stats["macs2 found in all"],
    stats["spp2 found in all"]]]) 

percent_tsv = "\n".join(["\t".join(k) for k in
    [sample_id,percent_stats["% spp-not-macs2"], percent_stats["% macs2-not-spp"], percent_stats["% macs2-not-spp2"],
    percent_stats["% spp2-not-macs2"], percent_stats["% spp-not-spp2"], percent_stats["% spp2-not-spp"],
    percent_stats["% spp_only"], percent_stats["% spp2_only"], percent_stats["% macs2_only"], percent_stats["% spp found in all"],
    percent_stats["% macs2 found in all"], percent_stats["% spp2 found in all"]]])

with open("bedtools_out/bedtools_stats.tsv", "w") as file:
    file.write(tsv)

with open("bedtools_out/bedtools_percent_stats.tsv", "w") as file:
    file.write(percent_tsv)

# Now create the heatmap or clustermap

fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(1, 1, 1)
fig.suptitle("Number of peaks by sample")
array = pd.DataFrame.from_csv('bedtools_out/bedtools_stats.tsv', sep='\t')
ax = sns.heatmap(array,square=True)
#cbar = plt.colorbar(heatmap)
#cbar.set_label("# of peaks", rotation=270)
fig.savefig(snakemake.output.heatmap)

fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(1, 1, 1)
fig.suptitle("Percentage of peaks distributed")
percent_array = pd.DataFrame.from_csv('bedtools_out/bedtools_percent_stats.tsv', sep='\t')
ax = sns.heatmap(percent_array,square=True)
#cbar = plt.colorbar(percent_heatmap)
#cbar.set_label("% of peaks", rotation=270)
fig.savefig(snakemake.output.percent_heatmap)


fig = plt.figure(figsize=(14,10))
clustermap = sns.clustermap(array)
#cbar = plt.colorbar(clustermap)
#cbar.set_label("# of peaks", rotation=270)
plt.savefig(snakemake.output.clustermap)

fig = plt.figure(figsize = (14,10))
percent_clustermap = sns.clustermap(percent_array)
plt.savefig(snakemake.output.percent_clustermap)
