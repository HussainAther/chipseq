import os
import yaml
from ago2 import helpers
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.spatial import distance
from scipy.cluster import hierarchy
from matplotlib import pyplot as plt
README = """
Clustered heatmaps of interval similarity using different algorithms.

Sort of surprising, but there is still no consensus on which algorithm to use
for assessing the similarity between sets of interval files (for example, to
see if peaks of one TF enriched for peaks of another).

I checked out the recent literature, and it seems that enrichment algorithms
seem to use either a permutation test (GAT, regionR, our previous method), an
existing statisitical test or metric (bedtools' fisher and jaccard tools), or
novel algorithms (GenometriCor and IntervalStats).

The algorithms
==============
I chose four different metrics, each of which calculate "enrichment" in quite
different ways.

GAT
---
To represent permutation testing I chose GAT [1, 2] because it is still
maintained (latest commit was in Sept), has great documentation and software
engineering, is written by one of the developers of pysam, and seems well
thought-out. It has functionality that our previous method doesn't, and seems
to have all of the features that regionR does. GAT reports a fold change
statistic along with an adjusted pval. I'm showing the log2 fold change value,
and setting anything that has an adjusted pval > 0.05 to zero. So the only
things with color are those that are significantly enriched/disenriched.  GAT
does calculate self-self values, but I'm not sure how to interpret them. The
docs don't have any suggestions. So here, I'm artificially setting self-self
comparisons to zero.

GAT also reports the fraction of nucleotides shared in the segments and
annotations ("query" and "reference" using IntervalStats nomenclature, below).
While this can be simply calculated using bedtools intersect, I'm using the GAT
output and showing it as a separate figure because it handles the domain
subsetting for free.

bedtools fisher
---------------
I'm using bedtools' fisher tool [3], which calcuates a Fisher's exact test on
a 2x2 contingency table of how many intervals overlap between two sets of
intervals and the intervals unique to each one.  The tricky part is figuring
out what to put in the lower-right corner of the FET table, which represents
"intervals not in either set". For that, you need to infer the total number of
intervals possible and therefore make assumptions about the average size and
number of intervals. See the docs at [3] for details, but basically the pvals
here are pretty inflated (lots of false positives) and as a result they
recommend doing permutation tests to follow up.

The fisher heatmap shows the -log10(pval) for the Fisher's exact test, which is
flipped to negative if the test indicates depletion rather than enrichment
(based on whether the "ratio" calculation is less than one or not).  Self-self
intersections and other cases where the two-tail pval is zero are set to the
global max value to avoid taking the log of zero.

bedtools jaccard
----------------
I'm also using bedtools' jaccard tool [4]. This is a simpler metric which is
the total bp in the intersection of two sets divided by the total bp in the
union of the two sets. That is, intersection divided by total among both sets.
The value varies from 0 (no overlap) to 1 (perfect overlap). I ended up
truncating the colorbar in the heatmap because there are many smaller values
that get swamped out by the perfectly-intersecting self-self intersections and
the promoter intersections. The Jaccard heatmap shows the Jaccard statistic for
each comparison.

(By the way: there is a bedtools reldist tool that I'm not using because
I don't know of a good way to summarize the resulting histograms into a single
value that can be put into a heatmap. It's inspired by GenometriCorr
R package).

Self-self intersections are defined by the Jaccard metric to be 1.0.

IntervalStats
-------------
Last, I'm using IntervalStats [5]. This is a more sophisticated algorithm that
actually gives a pval for the association of each interval in a query with
a reference. Unlike the other tools (with the exception of percent overlap),
*IntervalStats is asymmetric*. That is, the value for "AGO2 peaks vs LAD
borders" is different from the value for "LAD borders vs AGO2 peaks".  This
allows for more biologically-nuanced interpretation.

For example, consider a transcription factor's overlap with PolII where the TF
intersects a subset of promoter that also have PolII -- but there are many
places with PolII that don't have the TF.  With the TF as query and PolII as
reference, we'd get a strong overlap . But the reverse -- PolII as the query
and the TF as the reference -- would give a weaker association statistic
because . This makes biological sense: most TFs are found with PolII, but not
all PolII sites have the TF.

The values in the heatmap are "fraction of intervals with pval < 0.05", which
is what they recommend in the paper [5] as a summary statistic. IntervalStats
does calculate self-self intersections, but like GAT I'm not sure how to
interpret them so I'm setting them to 1.0.


The domain
==========
GAT calls it "workspace", IntervalStats calls it "domain": the region of the
genome where you're allowed to place an interval when calculating a statistic.
For now, I'm using uniquely mappable regions. I'm calculating these using the
GEM mapper [8,9], selecting those regions that are uniquely mappable with 50-bp
reads (Brian's lab uses GEM to calculate mappability as well).

Other options may be AGO2 peak regions, promoters, or LAD boundaries. For
example, using AGO2 peaks as a domain, we'd be asking the question "what is the
association within AGO2 peaks?".

Note that bedtools fisher and jaccard don't actually support
domains/workspaces; to sort of get around this I'm only keeping the intervals
that intersect the domain before feeding them to the respective algorithm. This
is kind of a hack but it's as close as I can get for now.

The plots
=========
For each algorithm I'm plotting a separate clustered heatmap. The exception is
GAT, which has two plots (one for log2foldchange, one for percent overlap; see
description above). I'm using correlation as the distance metric and the
linkage method is "average" (see refs 6 and 7). The values represented in the
heatmap are indicated in the colorbar label.

The interval sets
=================
Each set of intervals is converted to BED3 format and any features on
non-euchromatic chromosomes are removed.

* The AGO2-related tracks are 1-bp summits of chipseq peaks that I've expanded
  out to 500 bp on either side (so, 1-kb peaks)

* LAD means entire LAD.

* LAD-inner means the 3kb from LAD border towards the middle of the LAD, with
  a 100bp buffer around the border

* LAD-outer means the 3kb from LAD border away from the LAD, with a 100bp
  buffer around the border.

* siAGO2 and dsLamin up/down/unchanged use a threshold of padj < 0.1. They
  represent all unique isoform TSSs for each gene, extended by 1kb on either
  side and then merged so that adjacent promoter regions are fused together.
  The "promoters" track is constructed similarly

* insulators are from Wood 2011. The original sizes range from 250 to 2500 bp.
  To make them more consistent both internally and with the AGO2 peaks, I took
  the 1kb surrounding the midpoint of each.

* Phantom peaks are from [10]. Lee mentioned them to me the other day, so
  I converted them to dm6 and added them to the mix. They seem to cluster with
  insulators here, as they mention in that paper.

I've uploaded the intervals as a trackhub which you can check for debugging.
The URL to paste into the "Add hub" box is:

    http://helix.nih.gov/~dalerr/lei/interval-sets/interval-sets.hub.txt


Notes on interpretation
=======================

GAT and Fisher can show disenrichment. IntervalStats and Jaccard cannot.

IntervalStats is asymmetric. Read it as "association of QUERY with REFERENCE"
(query is row, reference is column).

IntervalStats is not very sensitive to size of intervals (which they show in
the paper). Other algorithms might or might not be sensitive. I tried to
mediate that by making the peak sizes uniform, but LADs still have varying
lengths.

Fisher doesn't seem that useful since it's detecting lots of highly significant
associations.

Jaccard seems to be clustering things that make sense. But there's no
assessment of significance; is "0.3" a high or low value? Everything needs to
be interpreted in biological context.

GAT has separate tools for directly comparing "segments" and "annotations",
which is what they call "query" and "reference". For example, it can address
the question, "does AGO2 associate with LAD borders more or less than CTCF?".
Think of it as comparing rows to each other to see which columns are
significantly different. I'm not doing this on everything since it's not clear
how best to summarize it, and the complexity of computation here will explode.
Keeping it simple for now.

I'm leaning towards using just GAT and IntervalStats.

Other technical notes
=====================
Notes for Ryan-of-the-future, stored here so they're all in one place.

* 23 files means 23 * 23 * 4 = 2116 independent runs (across GAT, IntervalStats,
 jaccard, fisher). This script handles the hierarchically-stored output data
 from those. Domain is also a dir level for when we need more tweaking.

* Data prepared with prepare_intervalstats_config.py. This in turn calls
 a lot of the caching machinery in the ago2/*_datastore modules. But since
 everything is to be compared together in these plots sometimes we have to do
 some post-processing (like adjusting peak widths; removing heterochromatic
 chroms).  prepare_intervalstats_config.py handles this and stores copies in
 prepare_intervalstats_data keyed by interval_config.yaml['beds'] key.

 NOTE: domains will have to be nested if running multiple domains on biowulf
 simultaneously!

* Metadata are also created by that script and stored in interval_config.yaml
  to be read by the snakefile.

* bigbeds are made too, and a trackhub is uploaded. See that script for
  details.

* Snakefile (plus slurm_cluster_config.yaml and WRAPPER_SLURM) can be used to
  submit jobs to biowulf using standard "sbatch --export=SNAKEFILE=Snakefile
  --time=$SLURM_TIME_LIMIT WRAPPER_SLURM". Careful of logfile explosion since
  I didn't set up the stderr/stdout merge for this.


References
==========
[1] GAT docs: http://gat.readthedocs.org/en/latest/
[2] GAT paper: http://bioinformatics.oxfordjournals.org/content/29/16/2046
[3] fisher docs: http://bedtools.readthedocs.org/en/latest/content/tools/fisher.html
[4] jaccard docs: http://bedtools.readthedocs.org/en/latest/content/tools/jaccard.html
[5] IntervalStats paper: http://bioinformatics.oxfordjournals.org/content/28/5/607.full
[6] graphical view of linkage methods: http://people.revoledu.com/kardi/tutorial/Clustering/Linkages.htm
[7] nice high-level interpretation of linkage methods: http://support.minitab.com/en-us/minitab/17/topic-library/modeling-statistics/multivariate/item-and-cluster-analyses/linkage-methods/
[8] GEM docs: http://algorithms.cnag.cat/wiki/The_GEM_library
[9] GEM paper: http://www.ncbi.nlm.nih.gov/pubmed/23103880
[10] Phantom peaks paper: http://www.ncbi.nlm.nih.gov/pubmed/26117547
"""
logger = helpers.get_logger()
figdir = helpers.figure_dir()

config = yaml.load(open('interval_config_s2.yaml'))

"""
IntervalStats columns:
1. Query interval
2. Closest reference interval
3. Length of query
4. Distance
5. Numerator
6. Denominator
7. p-value (quotient of above two)
"""

def dataframe_for_domain(domain, kind):
    df = []
    for query in config['beds']:
        logger.info(query)
        for reference in config['beds']:
            filename = os.path.join(
                config['output'],
                kind,
                domain,
                '{0}_vs_{1}.txt'.format(query, reference)
            )

            _df = pd.read_table(filename, comment='#')
            _df['query'] = query
            _df['reference'] = reference
            df.append(
                _df.ix[0].to_dict()
            )
    return pd.DataFrame(df)


METRIC = 'correlation'
METHOD = 'average'

conf = [
    ('GAT', 'l2fold'),
    ('IntervalStats', 'f_01'),
    ('fisher', 'two-tail'),
    ('jaccard', 'jaccard'),
    ('GAT', 'fractions'),
]


for kind, value in conf:
    df = dataframe_for_domain('uniquely-mappable', kind)
    vmin, vmax = None, None
    if kind == 'IntervalStats':
        piv = df.pivot(index='query', columns='reference', values='f_01')
        fill_piv = piv.fillna(-1)
        np.fill_diagonal(fill_piv.values, 1)
        units = 'fraction pvals < 0.05'
        title = 'IntervalStats'

    if kind == 'GAT' and value == 'l2fold':
        piv = df.pivot(index='query', columns='reference', values='l2fold')
        mask = df.pivot(index='query', columns='reference', values='qvalue')
        title = 'GAT foldchange'

        # anything not significant becomes l2fold = 0
        piv[mask>0.05] = 0
        fill_piv = piv
        np.fill_diagonal(fill_piv.values, 0)
        units = 'log2fold'

    if kind == 'GAT' and value == 'fractions':
        segment_frac = df.pivot(index='query', columns='reference', values='percent_overlap_size_track')
        annotation_frac = df.pivot(index='query', columns='reference', values='percent_overlap_size_annotation')
        mask = df.pivot(index='query', columns='reference', values='qvalue')
        piv = segment_frac
        lower_tri_mask = np.ones(piv.shape, dtype='bool')
        lower_tri_mask[np.tril_indices(len(piv))] = False
        piv[lower_tri_mask] = annotation_frac[lower_tri_mask]
        piv[mask > 0.05] = 0
        fill_piv = piv
        units = 'percentage overlap'
        title = 'GAT percentage nucleotide overlap'


    if kind == 'fisher':
        piv = df.pivot(index='query', columns='reference', values='two-tail')
        mask_left = df.pivot(index='query', columns='reference', values='left')
        mask_right = df.pivot(index='query', columns='reference', values='right')
        mask_ratio = df.pivot(index='query', columns='reference', values='ratio')

        flip = mask_ratio < 1
        piv = -np.log10(piv)
        piv[flip] *= -1
        mx = piv.replace([np.inf],0).max().max()
        mn = piv.replace([-np.inf],0).min().min()
        piv = piv.replace([np.inf], mx)
        piv = piv.replace([-np.inf], mn)
        fill_piv = piv.fillna(0)
        units = '-log10(pval)'
        title = 'Fisher'

    if kind == 'jaccard':
        piv = df.pivot(index='query', columns='reference', values='jaccard')
        fill_piv = piv
        units = 'Jaccard statistic'
        vmin, vmax = (0, .3)
        title = 'Jaccard'


    dist = distance.pdist(fill_piv.values, metric=METRIC)
    dist[np.isnan(dist)] = 0
    dist[dist < 0] = 0


    if METHOD == 'ward':
        vals = fill_piv
    else:
        vals = dist
    row_linkage = hierarchy.linkage(vals, method=METHOD)

    ind = hierarchy.dendrogram(row_linkage, no_plot=True)['leaves']
    a = sns.clustermap(fill_piv, row_linkage=row_linkage,
                       col_linkage=row_linkage, robust=True, linewidths=1,
                       vmin=vmin, vmax=vmax)
    for txt in a.ax_heatmap.get_xticklabels():
        txt.set_rotation(90)
    for txt in a.ax_heatmap.get_yticklabels():
        txt.set_rotation(0)

    a.cax.set_ylabel(units)

    fig = plt.gcf()
    fig.suptitle(title, weight='bold', size=20)
    fig.subplots_adjust(right=0.8, bottom=0.2)
    fn = os.path.join(figdir, kind + '_' + value + '.png')

    fig.savefig(fn)
    logger.info('saved %s' % fn)

helpers.write_readme(README)
