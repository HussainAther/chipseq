import matplotlib
import os
from matplotlib.mlab import bivariate_normal
matplotlib.use('Agg')
from snakemake import shell
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors
import scipy

from deeptools.heatmapper import heatmapper

def kind_from_filename(fn):
    """
    "profiles/reference-point/macs2_only/macs2_only_SRR567586.tsv" -> "macs2_only"
    """
    #return os.path.basename(os.path.dirname(fn))
    return str(fn).split("/")[2]

def sample_from_filename(fn, kind):
    """
    "matrices/reference-point/macs2_only/macs2_only_SRR567586.gz" -> "SRR567586"
    """
    return str(fn.split(".")[0].split("/")[3].replace(kind, "")[1:])

stats = {}
sample = ""
for filename in snakemake.input.profiles:
    kind = kind_from_filename(filename)
    sample = sample_from_filename(filename, kind)
    # if the file is empty, we need to decide what to do with it. For now, just
    # continue
    if os.stat(filename).st_size == 0:
        continue

    # use deeptools' parsing to handle loading the matrix. They actually
    # support some pretty complex features (concatentating upstream, gene body,
    # downstream; concatentating multiple sets of features). The snakefile is
    # running the simpler mode of a single set of features. So we can simply
    # grab the (0, 0) matrix.
    h = heatmapper()
    h.read_matrix_file(filename)

    # `matrix` is a numpy array, with one row per feature (in this context, one
    # row per peak) and one column for each bin
    matrix = h.matrix.get_matrix(0, 0)['matrix']

    # take the mean of each column
    y = matrix.mean(axis=0)

    # deeptools.Heatmapper.read_matrix_file also parses the parameters, which
    # is pretty nice. We can use that to build an x-axis.

    if snakemake.wildcards.scaling == 'reference-point':
        x = np.linspace(
            -h.parameters['upstream'],
            h.parameters['downstream'],
            matrix.shape[1]
        )
    elif snakemake.wildcards.scaling == 'scale-regions':
        x = np.arange(h.parameters['upstream'], h.parameters['body'], h.parameters['bin size'])
    else:
        raise ValueError('Unhandled scaling: {}'.format(snakemake.wildcards.scaling))

    # Since the values are all pretty different, it's hard to compare their
    # shapes on the same plot. So let's rescale them for plotting in
    # a different way.
    normed_y = y / y.max()

    peak_count = matrix.shape[0]

    # store those values in the stats dictionary for later
    stats[kind] = {
        'matrix': matrix,
        'x': x,
        'y': y,
        'normed_y': normed_y,
        'peak_count': peak_count,
    }

color_dict = {
    'macs2_only' : '#000000',
    'spp_only' : '#FF0000',
    'spp2_only' : '#808000',
    'spp-not-spp2' : '#00FF00',
    'spp2-not-spp' : '#008000',
    'macs2-not-spp' : '#00FFFF',
    'macs2-not-spp2' : '#0000FF',
    'spp-not-macs2' : '#800080',
    'spp2-not-macs2' : '#000080',
    'spp_found-in-all' : '#FF00FF',
    'macs2_found-in-all' : '#808080',
    'spp2_found-in-all' : '#800000'
}

fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)

for kind, d in stats.items():
    label = '{0} ({1})'.format(kind, d['peak_count'])
    ax1.plot(d['x'], d['y'], label=label, color=color_dict[kind])
    ax2.plot(d['x'], d['normed_y'], color=color_dict[kind])

ax1.set_ylabel('Average density')
ax2.set_ylabel('Normalized density')
ax1.set_xlabel('position')
ax2.set_xlabel('position')

handles, labels = ax1.get_legend_handles_labels()
fig.tight_layout()
fig.subplots_adjust(right=0.7, top=0.9)
if snakemake.wildcards.scaling == 'reference-point':
    fig.suptitle("Reference-point density profile of " + str(sample) + " from computeMatrix", fontsize='14')
elif snakemake.wildcards.scaling == 'scale-regions':
    fig.suptitle("Scale-regions density profile of " + str(sample) + " from computeMatrix", fontsize='14')
fig.legend(handles, labels, 'center right', fontsize='12')
fig.savefig(snakemake.output.lineplot)
