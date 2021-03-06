from matplotlib.colors import LogNorm
from matplotlib.mlab import bivariate_normal
from snakemake.shell import shell
import matplotlib
matplotlib.use('Agg')
import numpy
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import math
import pylab

# This file will look at the output 
# from the peak error rules and create output of 
# the false positive and false negative values and
# also of the matthews correlation coefficient values

def frange(x, y, jump):
    while x < y:
        yield x
        x+= jump

q_value = {} #keep track of the q_value stats dictionary for each q_value for each sample
q_list = []
for i in list(frange(10e-15, .8, ((.8-10e-15)/59))):
    if i != 1e-14:
        q_value[str(i)] = [0, 0, 0, 0, 0]
        q_list.append(str(i))

stats = {}
for i in q_value:
    stats[i] = q_value[i]

#print(stats)

def mcc(tp, tn, fp, fn): #calculate the matthews correlation coefficient
    num = tp*tn - fp*fn
    dem = math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    if dem == 0:
        dem = 1
    return num/dem

for file in snakemake.input.errors:
    f = open(file, "r")
    q_value = file.split("/")[3]
    for lineno, line in enumerate(f):
        if lineno>0:
            if line.split("\t")[11] == "correct":
                stats[q_value][0] += 1
            elif line.split("\t")[11] == "false negative":
                stats[q_value][1] += 1
            elif line.split("\t")[11] == "false positive":
                stats[q_value][3] += 1
            else:
                stats[q_value][2] +=1
    stats[q_value][4] += mcc(stats[q_value][0], stats[q_value][2], stats[q_value][3], stats[q_value][1])

fn = []
fp = []
fsum = []
tn = []
tp = []
tsum = []
mcc = []
for q_value in stats:
    fn.append(stats[q_value][1])
    fp.append(stats[q_value][3])
    fsum.append(stats[q_value][1] + stats[q_value][3])
    tn.append(stats[q_value][2])
    tp.append(stats[q_value][0])
    tsum.append(stats[q_value][2] + stats[q_value][0])
    mcc.append(stats[q_value][4])

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_title("False positives and false negatives of " + str(snakemake.params.sample))
ax.plot(fn, label="false negatives")
ax.plot(fp, label="false positives")
ax.plot(fsum, label="fp + fn")
ax.legend()
#plt.xticks(q_list)
ax.set_xlabel('q-value')
ax.set_ylabel('count')
fig.savefig(snakemake.output.fpfn)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_title("True positives and true negatives of " + str(snakemake.params.sample))
ax.plot(tn, label="true negatives")
ax.plot(tp, label="true positives")
ax.plot(tsum, label="tp + tn")
ax.legend()
#plt.xticks(q_list)
ax.set_xlabel('q-value')
ax.set_ylabel('count')
fig.savefig(snakemake.output.tptn)

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111)
ax.set_title("Matthews correlation coefficient of " + str(snakemake.params.sample))
ax.plot(mcc, label="mcc")
ax.legend()
#plt.xticks(q_list)
ax.set_xlabel('q-value')
ax.set_ylabel('count')
fig.savefig(snakemake.output.mcc)

