import matplotlib
import os
from matplotlib.colors import LogNorm
from matplotlib.mlab import bivariate_normal
matplotlib.use('Agg')
from snakemake import shell
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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
    return str(fn).split("/")[3].split(".")[0].replace(kind, "")[1:]

stats = {}

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
    stats[sample] = {
        'matrix': matrix,
        'x': x,
        'y': y,
        'normed_y': normed_y,
        'peak_count': peak_count,
    }

# Create the figures with random colors

fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)

for sample, d in stats.items():
    label = '{0} ({1})'.format(sample, d['peak_count'])
    ax1.plot(d['x'], d['y'], label=label)
    ax2.plot(d['x'], d['normed_y'])

ax1.set_ylabel('Average density')
ax2.set_ylabel('Normalized density')
ax1.set_xlabel('position')
ax2.set_xlabel('position')

handles, labels = ax1.get_legend_handles_labels()
fig.tight_layout()
fig.subplots_adjust(right=0.7, top=0.9)
if snakemake.wildcards.scaling == 'reference-point':
    fig.suptitle("Reference-point density profile of " + kind_from_filename(filename) + " from computeMatrix", fontsize='14')
elif snakemake.wildcards.scaling == 'scale-regions':
    fig.suptitle("Scale-regions density profile of " + kind_from_filename(filename) + " from computeMatrix", fontsize='14')
fig.legend(handles, labels, 'center right', fontsize='12')
fig.savefig(snakemake.output.lineplot)

# Now we'll create the figures but using only one color per cell-type with a legend that only has a single label for each color.

color_dict_celltype = {
    'el39_bg3_rump-5g4_R1' : '#FF0000',
    'el40_kc_rump-5g4_R1' : '#0000FF',
    'el41_bg3_rump-10c3_R1' : '#FF0000',
    'el42_kc_rump-10c3_R1' : '#0000FF',
    'lm06_kc_shep-gp_R1' : '#0000FF',
    'lm07_kc_shep-r5_R1' : '#0000FF',
    'lm08_kc_shep-r6_R1' : '#0000FF',
    'lm10_kc_shep-gp_R1' : '#0000FF',
    'lm32_kc_suhw_R1' : '#0000FF',
    'lm33_bg3_suhw_R1' : '#FF0000',
    'lm34_bg3_shep-gp_R1' : '#FF0000',
    'lm36_bg3_shep-r6_R1' : '#FF0000',
    'lm37_kc_mod-gp_R1' : '#0000FF',
    'lm38_kc_mod-rb_R1' : '#0000FF',
    'lm39_kc_rm62-rat_R1' : '#0000FF',
    'lm40_kc_rm62-rb_R1' : '#0000FF',
    'lm41_kc_rump-5g4_R1' : '#0000FF',
    'lm42_kc_rump-10c3_R1' : '#0000FF',
    'lm43_bg3_mod-gp_R1' : '#FF0000',
    'lm44_bg3_mod-rb_R1' : '#FF0000',
    'lm47_bg3_rump5g4_R1' : '#FF0000',
    'lm48_bg3_rump-10c3_R1' : '#FF0000',
    'sjl_cl8_cp190_R1' : '#008000',
    'sjl_cl8_mod_R1' : '#008000',
    'sjl_cl8_suhw_R1' : '#008000',
    'sjl_kc_cp190_R1' : '#0000FF',
    'sjl_mbn2_cp190_R1' : '#000000',
    'sjl_mbn2_mod_R1' : '#000000',
    'sjl_mbn2_shep_R1' : '#000000',
    'sjl_mbn2_suhw_R1' : '#000000',
    'sjl_s2_cp190_R1' : '#800080',
    'sjl_s2_mod_R1' : '#800080',
    'sjl_s2_shep_R1' : '#800080',
    'sjl_s3_cp190_R1' : '#800000',
    'sjl_s3_suhw_R1' : '#800000',
    'SRR567586' : '#FF0000',
    'dc92_bg3_ctcf-56' : '#000000',
    'dc94_bg3_ctcf-59' : '#000000'
    }

fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)

label_list = {} # Keep track of which cell type has been assigned a label for the legend. 
 # This way, the legend only shows one cell type for each color.

for sample, d in stats.items():
    if sample == "SRR567586":
        if "bg3" not in label_list:
            label_list["bg3"] = color_dict_celltype[sample]
            ax1.plot(d['x'], d['y'], label="bg3", color=color_dict_celltype[sample])
            ax2.plot(d['x'], d['normed_y'], color=color_dict_celltype[sample])
        else:
            ax1.plot(d['x'], d['y'], color=color_dict_celltype[sample])
            ax2.plot(d['x'], d['normed_y'], color=color_dict_celltype[sample])
    elif sample.split("_")[1] not in label_list:
        label_list[sample.split("_")[1]] = color_dict_celltype[sample]
        ax1.plot(d['x'], d['y'], label=sample.split("_")[1], color=color_dict_celltype[sample])
        ax2.plot(d['x'], d['normed_y'], color=color_dict_celltype[sample])
    else:
        ax1.plot(d['x'], d['y'], color=color_dict_celltype[sample])
        ax2.plot(d['x'], d['normed_y'], color=color_dict_celltype[sample])

ax1.set_ylabel('Average density')
ax2.set_ylabel('Normalized density')
ax1.set_xlabel('position')
ax2.set_xlabel('position')

handles, labels = ax1.get_legend_handles_labels()
fig.tight_layout()
fig.subplots_adjust(right=0.7, top=0.9)
if snakemake.wildcards.scaling == 'reference-point':
    fig.suptitle("Reference-point density profile of " + kind_from_filename(filename) + " from computeMatrix", fontsize='14')
elif snakemake.wildcards.scaling == 'scale-regions':
    fig.suptitle("Scale-regions density profile of " + kind_from_filename(filename) + " from computeMatrix", fontsize='14')
fig.legend(handles, labels, 'center right', fontsize='14')
fig.savefig(snakemake.output.lineplot_celltype)

# Now create the same lineplots but with a different color for each antibody. 
# When two antibodies have different numbers at the end, treat them as different antibodies. 

color_dict_antibody = {
    'el39_bg3_rump-5g4_R1' : '#e41a1c',
    'el40_kc_rump-5g4_R1' : '#e41a1c',
    'el41_bg3_rump-10c3_R1' : '#e41a1c',
    'el42_kc_rump-10c3_R1' : '#e41a1c',
    'lm06_kc_shep-gp_R1' : '#377eb8',
    'lm07_kc_shep-r5_R1' : '#377eb8',
    'lm08_kc_shep-r6_R1' : '#377eb8',
    'lm10_kc_shep-gp_R1' : '#377eb8',
    'lm32_kc_suhw_R1' : '#4daf4a',
    'lm33_bg3_suhw_R1' : '#4daf4a',
    'lm34_bg3_shep-gp_R1' : '#377eb8',
    'lm36_bg3_shep-r6_R1' : '#377eb8',
    'lm37_kc_mod-gp_R1' : '#984ea3',
    'lm38_kc_mod-rb_R1' : '#984ea3',
    'lm39_kc_rm62-rat_R1' : '#ff7f00',
    'lm40_kc_rm62-rb_R1' : '#ff7f00',
    'lm41_kc_rump-5g4_R1' : '#e41a1c',
    'lm42_kc_rump-10c3_R1' : '#e41a1c',
    'lm43_bg3_mod-gp_R1' : '#984ea3',
    'lm44_bg3_mod-rb_R1' : '#984ea3',
    'lm47_bg3_rump5g4_R1' : '#e41a1c',
    'lm48_bg3_rump-10c3_R1' : '#e41a1c',
    'sjl_cl8_cp190_R1' : '#ffff33',
    'sjl_cl8_mod_R1' : '#984ea3',
    'sjl_cl8_suhw_R1' : '#4daf4a',
    'sjl_kc_cp190_R1' : '#ffff33',
    'sjl_mbn2_cp190_R1' : '#ffff33',
    'sjl_mbn2_mod_R1' : '#984ea3',
    'sjl_mbn2_shep_R1' : '#377eb8',
    'sjl_mbn2_suhw_R1' : '#4daf4a',
    'sjl_s2_cp190_R1' : '#ffff33',
    'sjl_s2_mod_R1' : '#984ea3',
    'sjl_s2_shep_R1' : '#377eb8',
    'sjl_s3_cp190_R1' : '#ffff33',
    'sjl_s3_suhw_R1' : '#4daf4a',
    'SRR567586' : '#4daf4a',
    'dc92_bg3_ctcf-56' : '#f781bf',
    'dc94_bg3_ctcf-59' : '#f781bf'
    }

fig = plt.figure(figsize=(10, 8))
ax1 = fig.add_subplot(2, 1, 1)
ax2 = fig.add_subplot(2, 1, 2, sharex=ax1)

label_list = {} # Dictionary with antibodies as keys and the antibody color as sample

for sample, d in stats.items():
    if sample == "SRR567586":
        if "suhw" not in label_list:
            label_list["suhw"] = color_dict_antibody[sample]
            ax1.plot(d['x'], d['y'], label="suhw", color=color_dict_antibody[sample])
            ax2.plot(d['x'], d['normed_y'], color=color_dict_antibody[sample])
        else:
            ax1.plot(d['x'], d['y'], color=color_dict_antibody[sample])
            ax2.plot(d['x'], d['normed_y'], color=color_dict_antibody[sample])
    elif "rump" in sample:
        if "rump" not in label_list:
            label_list["rump"] = color_dict_antibody[sample]
            ax1.plot(d['x'], d['y'], label="rump", color=color_dict_antibody[sample])
            ax2.plot(d['x'], d['normed_y'], color=color_dict_antibody[sample])
        else:
            ax1.plot(d['x'], d['y'], color=color_dict_antibody[sample])
            ax2.plot(d['x'], d['normed_y'], color=color_dict_antibody[sample])
    elif sample.split("_")[2].split("-")[0] not in label_list:
        label_list[sample.split("_")[2].split("-")[0]] = color_dict_antibody[sample]
        ax1.plot(d['x'], d['y'], label=sample.split("_")[2].split("-")[0], color=color_dict_antibody[sample])
        ax2.plot(d['x'], d['normed_y'], color=color_dict_antibody[sample])
    else:
        ax1.plot(d['x'], d['y'], color=color_dict_antibody[sample])
        ax2.plot(d['x'], d['normed_y'], color=color_dict_antibody[sample])

ax1.set_ylabel('Average density')
ax2.set_ylabel('Normalized density')
ax1.set_xlabel('position')
ax2.set_xlabel('position')

handles, labels = ax1.get_legend_handles_labels()
fig.tight_layout()
fig.subplots_adjust(right=0.7, top=0.9)
if snakemake.wildcards.scaling == 'reference-point':
    fig.suptitle("Reference-point density profile of " + kind_from_filename(filename) + " from computeMatrix", fontsize='14')
elif snakemake.wildcards.scaling == 'scale-regions':
    fig.suptitle("Scale-regions density profile of " + kind_from_filename(filename) + " from computeMatrix", fontsize='14')
fig.legend(handles, labels, 'center right', fontsize='12')
fig.savefig(snakemake.output.lineplot_antibody)
