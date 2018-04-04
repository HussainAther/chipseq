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

stats = {} # Keep track of the jaccard index for each file

# Iterate through each input file and collect the data

for file in snakemake.input:
    f = open(file, "r")
    for lineno, line in enumerate(f):
        if lineno == 1:
            stats[file] = line.split("\t")[2]

# Now create the tsv with those jaccard indices

file = snakemake.input

tsv = "\n".join(["\t".join(k) for k in 
    [file, stats[file]]])

with open(snakemake.output.stats_tsv, "w") as file:
    file.write(tsv)
