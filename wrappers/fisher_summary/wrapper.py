from snakemake import shell
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

# Summarize and visualize the Fisher's exact test results for each sample

stats = {} # Keep track of the 2x2 grid for each Fisher's exact test

for file in snakemake.input:
    stats[file] = [0, 0, 0, 0] #in a and b, in a not b, in b not a, not in a or b 

pvalue = {} # Keep track of the p_value for Fisher's exact test

for file in snakemake.input:
    pvalue[file] = [0, 0, 0, 0] # left, right, two-tail, ratio

# Iterate through each input file and collect the data

for file in snakemake.input:
    f = open(file, "r")
    for lineno, line in enumerate(f):
        if lineno == 8:
            stats[file][0] = line.split("|")[1].replace(" ","") 
            stats[file][1] = line.split("|")[2].replace(" ","") 
        if lineno == 9:
            stats[file][2] = line.split("|")[1].replace(" ","") 
            stats[file][3] = line.split("|")[2].replace(" ","") 
        if lineno == 13:
            pvalue[file][0] = line.split("\t")[0]
            pvalue[file][1] = line.split("\t")[1]
            pvalue[file][2] = line.split("\t")[2]
            pvalue[file][3] = line.split("\t")[3]

# Now create the tsv files with these statistics

file = snakeamke.input

stats_tsv = "\n".join(["\t".join(k) for k in 
    [file, stats[file][0], stats[file][1], stats[file][2], stats[file][3]]])

pvalue_tsv = "\n".join(["\t".join(k) for k in 
    [file, p_value[file][0], p_value[file][1], p_value[file][2], p_value[file][3]]])

with open(snakemake.output.stats_tsv, "w") as file:
    file.write(stats_tsv)

with open(snakemake.output.pvalue_tsv, "w") as file:
    file.write(pvalue_tsv)


# Now create the heatmap or clustermap

plt.title("Fisher's exact test stats")
plt.xlabel("sample")
plt.ylabel("counts")
array = pd.DataFrame.from_csv('bedtools_out/fisher/fisher_stats.tsv', sep='\t')
fig = plt.figure(figsize=(12,5))
fig = sns.heatmap(array,square=True)
#cbar = plt.colorbar(heatmap)
#cbar.set_label("# of peaks", rotation=270)
plt.savefig(snakemake.output.stats_summary)

plt.title("Fisher's exact test p_values")
plt.xlabel("sample")
plt.ylabel("value")
percent_array = pd.DataFrame.from_csv('bedtools_out/fisher/fisher_pvalue.tsv', sep='\t')
fig = plt.figure(figsize=(12,5))
fig = sns.heatmap(percent_array,square=True)
#cbar = plt.colorbar(percent_heatmap)
#cbar.set_label("% of peaks", rotation=270)
plt.savefig(snakemake.output.pvalue_summary)
