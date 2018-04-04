import matplotlib
matplotlib.use('Agg')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()
columns=[1,4]

spp_array=pd.read_csv('peak_out/spp/SRR567586/SRR567586.bed.narrowPeak.sorted', sep="\t", usecols=columns)
macs2_array=pd.read_csv('peak_out/macs2/SRR567586/SRR567586_peaks.narrowPeak.sorted', sep="\t",usecols=columns)
print(spp_array)

frames = [spp_array, macs2_array]
result = pd.concat(frames)

fig=plt.figure()
result_heatmap=sns.heatmap(spp_array)
fig.savefig('bedline_out.png')
