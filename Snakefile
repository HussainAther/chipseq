from snakemake import shell
import itertools
import pandas as pd
import os
import csv

metadata = pd.read_csv('config/metadata.tsv', header=0, sep='\t')
metadata = metadata+"_R1"
temp_dict = pd.DataFrame.to_dict(metadata, orient='list')

d = dict(zip(temp_dict['ip'], temp_dict['input'])) # dictionary with ip as keys and input as values
peakcallers = ['macs2','spp', 'spp2']
d["SRR567586"] = "SRR567585"
d["dc92_bg3_ctcf-56"] = "dc91_bg3_input"
d["dc94_bg3_ctcf-59"] = "dc93_bg3_input"

d_rev = {} # dictionary with input as keys and listed ip as values
for i in d:
    if d[i] not in d_rev:
        d_rev[d[i]] = [str(i)]
    else:
        d_rev[d[i]].append(i)

sample = [str(i) for i in d.keys()]
inp_sample = [str(i) for i in d_rev.keys()]
all_sample = sample + inp_sample
#sample = ['lm32_kc_suhw_R1', 'sjl_cl8_suhw_R1']

peakcallers = ['macs2', 'spp', 'spp2']
unique_to_one_peakcaller = [pc + '_only' for pc in peakcallers]
shared_across_all_peakcallers = [pc + '_found-in-all' for pc in peakcallers]
#pairwise_unique = ['{0}-not-{1}'.format(i, j) for i, j in itertools.permutations(peakcallers, 2)]
kinds = unique_to_one_peakcaller + shared_across_all_peakcallers

chr = ["chr2L", "chr2R", "chr3L", "chr3R", "chr4"]

kinds = [
   'macs2_only',
   'spp_only',
   'spp2_only',
   'macs2_found-in-all',
   'spp_found-in-all',
   'spp2_found-in-all',
   'macs2-not-spp',
   'macs2-not-spp2',
   'spp-not-macs2',
   'spp-not-spp2',
   'spp2-not-macs2',
   'spp2-not-spp'
]


# this function lets you return a range in python even with floats
def frange(x, y, jump):
    while x < y:
        yield x
        x += jump

q_value = []
for i in list(frange(10e-15,.8,((.8-10e-15)/59))):
    if i != 1e-14:
        q_value.append(i)
q_value = list(q_value)

targets = expand(
    '{plot_type}/{subprogram}/{kind}/{kind}_{sample}.pdf',
    plot_type=['heatmaps', 'profiles'],
    kind=kinds,
    subprogram=['scale-regions', 'reference-point'],
    sample=sample)


peaks = (
    expand(
        "peak_out/spp/{ip}/{ip}.bed",
        ip = sample,
    ) +
    expand(
        "peak_out/spp2/{ip}/{ip}.bed",
        ip = sample,
    ) +
    expand(
        "peak_out/{peakcallers}/{sample}/{sample}.bed.narrowPeak.sorted",
        sample = sample,
        peakcallers = ['spp', 'spp2']
        ) +
    expand(
        "peak_out/macs2/{sample}/{sample}_peaks.narrowPeak.sorted",
        sample = sample
    )
)


deeptools = (
    expand(
        "bigwigcompare/{ip}.bedgraph",
        ip = sample,
    ) +
    expand(
        "bigwigcompare/{ip}_ratio.bedgraph",
        ip = sample,
    ) +
    expand(
        "fingerprints/{ip}.pdf",
        ip = sample,
    )
)

fastqc = (
    expand(
        "fastqc/{ip}_fastqc.zip",
        ip = sample,
    ) +
    expand(
        "fastqc/{inp}_fastqc.zip",
        inp = inp_sample,
    ) +
    expand(
        "multiqc/{ip}",
        ip = sample,
    ) +
    expand(
        "multiqc/{inp}",
        inp = inp_sample,
    ) +
    ["multiqc/all"]
)

reports = (
#    expand(
#        "reports/st.html",
#        ip = sample,
#    ) +
#    expand(
#        "reports/exp.html",
#        ip = sample,
#    ) +
#    ["config/metadata.rds"] +
    [expand("profiles/reference-point/{kind}/{kind}_{sample}.tsv", kind = kinds, sample = sample),
    expand("profiles/scale-regions/{kind}/{kind}_{sample}.tsv", kind = kinds, sample = sample)] +
    [   "bedtools_out/heatmap.png",
        "bedtools_out/percent_heatmap.png",
        "bedtools_out/clustermap.png",
        "bedtools_out/percent_clustermap.png"] +
   [ expand("bedtools_out/lineplots/reference-point/lineplot_by_kind_{kind}.png", kind = kinds),
    expand("bedtools_out/lineplots/scale-regions/lineplot_by_kind_{kind}.png", kind = kinds),
    expand("bedtools_out/lineplots/reference-point/lineplot_by_kind_cell_{kind}.png", kind = kinds),
    expand("bedtools_out/lineplots/scale-regions/lineplot_by_kind_cell_{kind}.png", kind = kinds),
    expand("bedtools_out/lineplots/reference-point/lineplot_by_kind_antibody_{kind}.png", kind = kinds),
    expand("bedtools_out/lineplots/scale-regions/lineplot_by_kind_antibody_{kind}.png", kind = kinds)] +
    expand("bedtools_out/lineplots/reference-point/lineplot_by_sample_{sample}.png", sample = sample) +
    expand("bedtools_out/lineplots/scale-regions/lineplot_by_sample_{sample}.png", sample = sample)
)

corr = (
    expand("correlations/heatmap_{inp}_spp.png",
        inp = sample,
    ) +
    expand("correlations/scat_{inp}_spp.png",
        inp = sample,
    ) +
    expand("correlations/heatmap_{inp}_macs2.png",
        inp = sample,
    ) +
    expand("correlations/scat_{inp}_macs2.png",
        inp = sample,
    )
)


colocalization = (
    expand(
         "bedtools_out/jaccard/jaccard_{sample}.txt",
            sample=sample
    ) +
    expand(
        "bedtools_out/fisher/fisher_{sample}.txt",
            sample=sample
    ) +
#    ["bedtools_out/jaccard/jaccard_summary.png"] +
#    ["bedtools_out/fisher/fisher_stats_summary.png", 
#        "bedtools_out/fisher/fisher_pvalue_summary.png"] +
    expand(
        "gat/{sample}_seg_spp",
        sample=sample,
    ) +
    expand(
        "gat/{sample}_seg_macs2",
        sample=sample,
    ) +
    [ "colocalization/colocalization_heatmap.png",
      "gat/gat_heatmap.png",
      "intervalstats/intervalstats_heatmap.png"
    ]

)

ML = (
    expand(
        "peakerror/macs2/{sample}/{q_value}/{chr}/{sample}_errors.tsv",
        q_value = q_value,
        sample = sample,
        chr = chr
    ) +
    expand(
        "peakerror/macs2/{sample}/{q_value}/{chr}/{sample}_summary.tsv",
        q_value = q_value,
        sample = sample,
        chr = chr
    ) +
    [expand("peakerror/macs2/{sample}/mcc.png", sample = sample) +
    expand("peakerror/macs2/{sample}/fpfpn.png", sample = sample) +
    expand( "peakerror/macs2/{sample}/tptn.png", sample = sample)] +
    expand(
        "peak_out/macs2/{ip}/{q_value}/{ip}_peaks.broadPeak",
        ip = sample,
        q_value = q_value)
)

localrules: regex

rule all:
    input: peaks + deeptools + fastqc + reports + corr + targets 

rule fasta_fix:
    input: "data/dm6.fa"
    output: "data/dm6.fixed.fa"
    script: "tools/fasta_parser.py"

rule bowtie2_index:
    input:
        "data/dm6.fixed.fa"
    output:
        "data/dm6.fixed.1.bt2",
        "data/dm6.fixed.rev.1.bt2"
    params:
        prefix = "mapped_reads/dm6.fixed"
    wrapper:
        'file:wrappers/bowtie2_index'

rule bowtie2_map:
    input:
        fastq = "data/fastqs/{sample}.fastq.gz",
        index = "mapped_reads/dm6.fixed.1.bt2"
    output:
        "mapped_reads/{sample}.bam",
    params:
        prefix = "mapped_reads/dm6.fixed",
    log:
        "logs/bowtie/{sample}.log"
    wrapper:
        'file:wrappers/bowtie2_map'

rule fastqc:
    input:
        "data/fastqs/{sample}.fastq.gz"
    output:
        "fastqc/{sample}_fastqc.zip"
    wrapper:
        "file:wrappers/fastqc"

rule multiqc:
    input:
        fastqc = "fastqc/{sample}_fastqc.zip",
        markduplicates = "dedup/{sample}.bam"
    output:
        "multiqc/{sample}"
    wrapper:
        "file:wrappers/multiqc"

rule multiqc_all:
    input:
        fastqc = expand("fastqc/{sample}_fastqc.zip", sample=sample),
        markduplicates = expand("dedup/{sample}.bam", sample=sample)
    output:
        "multiqc/all"
    wrapper:
        "file:wrappers/multiqc_all"

rule bam_trunc:
    input:
        "dedup/{sample}.bam"
    output:
        "dedup/{sample}.bam.tr"
    wrapper:
        "file:wrappers/bam_trunc"

rule chipqc1:
    input:
        reads =  "dedup/{sample}.bam.tr",
        spp_peaks = "peak_out/spp/{sample}/{sample}.bed.narrowPeak.sorted.tr",
        macs2_peaks = "peak_out/macs2/{sample}/{sample}_peaks.narrowPeak.sorted.tr",
        spp2_peaks = "peak_out/spp2/{sample}/{sample}.bed.narrowPeak.sorted.tr"
    output:
       spp = "dedup/rds/spp_{sample}.rds",
       macs2 = "dedup/rds/macs2_{sample}.rds",
       spp2 = "dedup/rds/spp2_{sample}.rds"
    wrapper:
        "file:wrappers/chipqc1"

rule chipqc2:
    input:
        metadata = "config/metadata.tsv"
    output:
        metadata = "config/metadata.rds",
        st_html = "reports/st.html",
        exp_html = "reports/exp.html"
    wrapper:
        "file:wrappers/chipqc2"

rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam",
    output:
        "sorted_reads/{sample}.bam",
    shell:
        "samtools sort {input} -o {output} -O BAM"

rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index -b {input} {output}"

rule macs2: #run macs2 with unsupervised training
    input:
        ip = "dedup/{sample}.bam",
        inp = lambda wc: "dedup/{0}.bam".format(d[wc.sample])
    output:
        narrowPeak = "peak_out/macs2/{sample}/{sample}_peaks.narrowPeak",
        #bed = "peak_out/macs2/{sample}/{sample}.bed"
    params:
        sample_prefix = "peak_out/macs2/{sample}"
    wrapper:
        'file:wrappers/macs2'

# rule macs2_v2: #used in the supervised training process to hyperparameterize the macs options
#     input:
#         ip = "dedup/{sample}.bam",
#         inp = lambda wc: "dedup/{0}.bam".format(d[wc.sample])
#     output:
#         broadPeak = "peak_out/macs2/{sample}/{q_value}/{sample}_peaks.broadPeak",
#     params:
#         qval_prefix = "peak_out/macs2/{sample}/{q_value}",
#         sample_prefix = "peak_out/macs2/{sample}/{q_value}/{sample}",
#         qval = "{q_value}",
#     wrapper:
#         'file:wrappers/macs2_v2'

rule spp:
    input:
        ip = 'dedup/{sample}.bam',
        inp = lambda wc: "dedup/{0}.bam".format(d[wc.sample])
    output:
        bed = 'peak_out/spp/{sample}/{sample}.bed',
        narrowPeak = "peak_out/spp/{sample}/{sample}.bed.narrowPeak"
    params:
        prefix="peak_out/spp/{sample}",
        tecfilter = "TRUE"
    wrapper:
        'file:wrappers/spp'

rule spp2:
    input:
        ip = "dedup/{sample}.bam",
        inp = lambda wc: "dedup/{0}.bam".format(d[wc.sample])
    output:
        bed = "peak_out/spp2/{sample}/{sample}.bed",
        narrowPeak = "peak_out/spp2/{sample}/{sample}.bed.narrowPeak"
    params:
        prefix = "peak_out/spp2/{sample}",
        tecfilter = "FALSE"
    wrapper:
        "file:wrappers/spp"

rule ccat:
    input:
        ip = "dedup/{sample}.bam",
        inp = lambda wc: "dedup/{0}.bam".format(d[wc.sample])
    output:
        "peak_out/ccat/{sample}/{sample}.signficant.peak",
        "peak_out/ccat/{sample}/{sample}.significant.region"
    params:
        size = "data/dm6.chr.size",
        config = "wrappers/ccat/config.yml",
        project_name = "peak_out/ccat/{sample}"
    wrapper:
        "file:wrappers/ccat"

rule peaksegjoint:
    input:
        ip = "dedup/{sample}.bam",
        inp = lambda wc: "dedup/{0}.bam".format(d[wc.sample])
    output:
        bed = "peak_out/peaksegjoint/{sample}/{sample}.bed",
    wrapper:
        "file:wrappers/peaksegjoint"

rule regex:
    input:
        "peak_out/macs2/{sample}/{q_value}/{sample}_peaks.broadPeak"
    output:
        "peak_out/macs2/{sample}/{q_value}/{chr}/{sample}_peaks.broadPeak.threecolumns"
    params:
        chr = "{chr}"
    shell:
        """awk '{{OFS="\t"; print $1,$2,$3}}' {input} | grep '^{params.chr}.*$' -u | sort -k2 -n | awk '!seen[$0]++' > {output}"""

rule peakerror_compute:
    input:
        "peak_out/macs2/{sample}/{q_value}/{chr}/{sample}_peaks.broadPeak.threecolumns",
    output:
        "peakerror/macs2/{sample}/{q_value}/{chr}/{sample}_errors.tsv"
    params:
        labels = "peakerror/{chr}labels.bed",
    shell:
        "Rscript peakerror/PeakError_compute.R {input} {params.labels} > {output}"

rule peakerror_summarize:
    input:
        "peakerror/macs2/{sample}/{q_value}/{chr}/{sample}_errors.tsv"
    output:
        "peakerror/macs2/{sample}/{q_value}/{chr}/{sample}_summary.tsv"
    shell:
        "Rscript peakerror/PeakError_summarize.R {input} > {output}"

rule peakerror_visualize:
    input:
        errors = expand("peakerror/macs2/{{sample}}/{q_value}/{chr}/{{sample}}_errors.tsv",  q_value=q_value, chr=chr),
        summary = expand("peakerror/macs2/{{sample}}/{q_value}/{chr}/{{sample}}_summary.tsv", q_value=q_value,chr=chr)
    output:
        mcc = "peakerror/macs2/{sample}/mcc.png",
        fpfn = "peakerror/macs2/{sample}/fpfpn.png",
        tptn = "peakerror/macs2/{sample}/tptn.png"
    params:
        sample = "{sample}"
    wrapper:
        "file:wrappers/peakerror_visualize"

rule pepr:
    input:
        ip = 'dedup/{sample}.bam',
        inp = lambda wc: "dedup/{0}.bam".format(d[wc.sample])
    output:
        'peak_out/pepr/{sample}/{sample}.bed'
    wrapper:
        'file:wrappers/pepr'


rule mark_duplicates:
    input:
        "sorted_reads/{sample}.bam",
    output:
        bam = "dedup/{sample}.bam",
        metrics = "dedup/metrics/{sample}.metrics.txt",
    log:
        "logs/picard/dedup/{sample}.log",
    params:
        "REMOVE_DUPLICATES=true"
    wrapper:
        "file:wrappers/mark_duplicates"

rule bedtools_all:
    input:
        beds = ["peak_out/spp/{sample}/{sample}.bed.narrowPeak.sorted"] +
        ["peak_out/macs2/{sample}/{sample}_peaks.narrowPeak.sorted"] +
        ["peak_out/spp2/{sample}/{sample}.bed.narrowPeak.sorted"]
    output:
        'bedtools_out/{kind}/{kind}_{sample}.txt'
    params:
        kind = "{kind}"
    script:
        "file:wrappers/bedtools_all/wrapper.py"

rule sort:
    input:
        beds = "peak_out/{peak_caller}/{sample}/{sample}{ext}.narrowPeak"
    output:
        sorted = "peak_out/{peak_caller}/{sample}/{sample}{ext}.narrowPeak.sorted"
    shell:
        "sort -k 1,1 -k2,2n {input.beds} > {output.sorted}"

rule bedtools_truncate:
    input:
        spp_sorted = "peak_out/spp/{sample}/{sample}.bed.narrowPeak.sorted",
        macs2_sorted = "peak_out/macs2/{sample}/{sample}_peaks.narrowPeak.sorted",
        spp2_sorted = "peak_out/spp2/{sample}/{sample}.bed.narrowPeak.sorted",
        genome = "data/dm6.fixed.fa"
    output:
        spp_sorted_tr = "peak_out/spp/{sample}/{sample}.bed.narrowPeak.sorted.tr",
        macs2_sorted_tr = "peak_out/macs2/{sample}/{sample}_peaks.narrowPeak.sorted.tr",
        spp2_sorted_tr = "peak_out/spp2/{sample}/{sample}.bed.narrowPeak.sorted.tr"
    wrapper:
        "file:wrappers/bedtools_truncate"

rule map:
    input:
        spp_sorted = expand("peak_out/spp/{sample}/{sample}.bed.narrowPeak.sorted.tr", sample=d.keys()),
        macs2_sorted = expand("peak_out/macs2/{sample}/{sample}_peaks.narrowPeak.sorted.tr", sample=d.keys()),
        spp2_sorted = expand("peak_out/spp2/{sample}/{sample}.bed.narrowPeak.sorted.tr", sample=d.keys()),
    output:
        heatmap = "bedtools_out/heatmap.png",
        percent_heatmap = "bedtools_out/percent_heatmap.png",
        clustermap = "bedtools_out/clustermap.png",
        percent_clustermap = "bedtools_out/percent_clustermap.png",
    script:
        "file:wrappers/map/wrapper.py"


rule lineplots_by_kind:
    input:
        profiles = expand("matrices/{{scaling}}/{{kind}}/{{kind}}_{sample}.gz", sample=sample),
    output:
        lineplot = "bedtools_out/lineplots/{scaling}/lineplot_by_kind_{kind}.png",
        lineplot_celltype = "bedtools_out/lineplots/{scaling}/lineplot_by_kind_cell_{kind}.png",
        lineplot_antibody = "bedtools_out/lineplots/{scaling}/lineplot_by_kind_antibody_{kind}.png"
    script:
        "file:wrappers/lineplots_by_kind/wrapper.py"

rule lineplots_by_sample:
    input:
        profiles = expand("matrices/{{scaling}}/{kind}/{kind}_{{sample}}.gz", kind=kinds),
    output:
        lineplot = "bedtools_out/lineplots/{scaling}/lineplot_by_sample_{sample}.png",
    script:
        "file:wrappers/lineplots_by_sample/wrapper.py"

rule jaccard:
    input:
        spp_sorted = "peak_out/spp/{sample}/{sample}.bed.narrowPeak.sorted.tr",
        macs2_sorted = "peak_out/macs2/{sample}/{sample}_peaks.narrowPeak.sorted.tr",
    output:
        jaccard = "bedtools_out/jaccard/jaccard_{sample}.txt",
    wrapper:
        "file:wrappers/jaccard"

rule jaccard_summary:
    input: 
        "bedtools_out/jaccard/jaccard_{sample}.txt"
    output: 
        summary = "bedtools_out/jaccard/jaccard_summary.png",
        stats_tsv = "bedtools_out/jaccard/jaccard_stats.tsv",
    wrapper: 
        "file:wrappers/jaccard_summary"

rule fisher:
    input:
        spp_sorted = "peak_out/spp/{sample}/{sample}.bed.narrowPeak.sorted.tr",
        macs2_sorted = "peak_out/macs2/{sample}/{sample}_peaks.narrowPeak.sorted.tr",
    output:
        fisher = "bedtools_out/fisher/fisher_{sample}.txt",
    wrapper:
        "file:wrappers/fisher"

rule fisher_summary:
    input:
        expand("bedtools_out/fisher/fisher_{sample}.txt", sample=sample)
    output:
        sgats_summary = "bedtools_out/fisher/fisher_stats_summary.png",
        pvalue_summary = "bedtools_out/fisher/fisher_pvalue_summary.png",
        stats_tsv = "bedtools_out/fisher/fisher_stats.tsv",
        p_value_tsv = "bedtools_out/fisher/fisher_pvalue.tsv"
    wrapper:
        "file:wrappers/fisher_summary"

rule bamcoverage:
    input:
        "dedup/{sample}.bam",
    output:
        "coverage/{sample}_coverage.bw"
    wrapper:
        "file:wrappers/bamcoverage"

rule bigwigcompare:
    input:
        ip = "coverage/{sample}_coverage.bw",
        inp = lambda wc:"coverage/{0}_coverage.bw".format(d[wc.sample]),
    output:
        log2 = "bigwigcompare/{sample}.bedgraph",
        ratio = "bigwigcompare/{sample}_ratio.bedgraph"
    wrapper:
        "file:wrappers/bigwigcompare"

#rule matrix_plot1:
#    input:
#        ip_bam = "dedup/{sample}.bam",
#        inp_bam = lambda wc: "dedup/{0}.bam".format(d[wc.sample]),
#        ip = "coverage/{sample}_coverage.bw",
#        inp = lambda wc: "coverage/{0}_coverage.bw".format(d[wc.sample]),
#        spp_region = "peak_out/spp/{sample}/{sample}.bed.narrowPeak.sorted.tr",
#        spp2_region = "peak_out/spp2/{sample}/{sample}.bed.narrowPeak.sorted.tr",
#        macs2_region = "peak_out/macs2/{sample}/{sample}_peaks.narrowPeak.sorted.tr",
# #      pepr_region="peak_out/pepr/{sample}/{sample}.bed"
#        spp_only = "bedtools_out/spp_only/spp_only_{sample}_spp_macs2.txt",
#        macs2_only = "bedtools_out/macs2_only/macs2_only_{sample}_spp_macs2.txt",
#        spp2_only = "bedtools_out/spp2_only/spp2_only_{sample}_spp2_macs2.txt",
#        overlap_spp = "bedtools_out/overlap/overlap_{sample}_spp_macs2.txt",
#        overlap_spp2 = "bedtools_out/overlap/overlap_{sample}_spp2_macs2.txt"
#    output:
#        fingerprints = "fingerprints/{sample}.pdf",
#    params:
#        sample = "{sample}",
#        sub = "reference-point"
#    wrapper:
#        "file:wrappers/matrix_plot1"

rule compute_matrix:
    input:
        ip = "coverage/{sample}_coverage.bw",
        inp = lambda wc: "coverage/{0}_coverage.bw".format(d[wc.sample]),
        regions = 'bedtools_out/{kind}/{kind}_{sample}.txt',
    output:
        gzipped_matrix = 'matrices/{subprogram}/{kind}/{kind}_{sample}.gz',
        tsv_out = 'matrices/{subprogram}/{kind}/{kind}_{sample}.tsv',
    conda: 'wrappers/compute_matrix/environment.yaml'
    params:
        subprogram = "{subprogram}"
    script:
        "file:wrappers/compute_matrix/wrapper.py"

rule plot_heatmap:
    input:
        gzipped_matrix = 'matrices/{subprogram}/{kind}/{kind}_{sample}.gz'
    output:
        pdf='heatmaps/{subprogram}/{kind}/{kind}_{sample}.pdf'
    conda: 'wrappers/plot_heatmap/environment.yaml'
    params:
        kind = "{kind}",
        sample = "{sample}",
        subprogram = "{subprogram}"
    script:
        "file:wrappers/plot_heatmap/wrapper.py"

rule plot_profile:
    input:
        gzipped_matrix = 'matrices/{subprogram}/{kind}/{kind}_{sample}.gz'
    output:
        pdf = "profiles/{subprogram}/{kind}/{kind}_{sample}.pdf",
        tsv = "profiles/{subprogram}/{kind}/{kind}_{sample}.tsv",
    conda: 'wrappers/plot_profile/environment.yaml'
    script:
        "file:wrappers/plot_profile/wrapper.py"

rule plot_fingerprint:
    input:
        ip_bam = "dedup/{sample}.bam",
        inp_bam = lambda wc: "dedup/{0}.bam".format(d[wc.sample]),
    output:
        "fingerprints/{sample}.pdf"
    conda: "wrappers/matrix_plot1/environment.yaml"
    shell:
        "plotFingerprint -b {input.ip_bam} "
        "{input.inp_bam} "
        "-plot {output}"

rule correlations:
    input:
        ip = "coverage/{sample}_coverage.bw",
        inp = lambda wc: "coverage/{0}_coverage.bw".format(d[wc.sample]),
        spp_region = "peak_out/spp/{sample}/{sample}.bed.narrowPeak.sorted.tr",
        spp2_region = "peak_out/spp2/{sample}/{sample}.bed.narrowPeak.sorted.tr",
        macs2_region = "peak_out/macs2/{sample}/{sample}_peaks.narrowPeak.sorted.tr",
    output:
        mbs_spp = "multiBigwigSummary/{sample}_spp",
        mbs_macs2 = "multiBigwigSummary/{sample}_macs2",
        correlations_heat_spp = "correlations/heatmap_{sample}_spp.png",
        correlations_scat_spp = "correlations/scat_{sample}_spp.png",
        correlations_heat_macs2 = "correlations/heatmap_{sample}_macs2.png",
        correlations_scat_macs2 = "correlations/scat_{sample}_macs2.png",
    wrapper:
        "file:wrappers/correlations"

rule gat:
    input:
        spp = "peak_out/spp/{sample}/{sample}.bed.narrowPeak.sorted.tr",
        macs2 = "peak_out/macs2/{sample}/{sample}_peaks.narrowPeak.sorted.tr",
        spp2 = "peak_out/spp2/{sample}/{sample}.bed.narrowPeak.sorted.tr"
    output:
        spp_nuc = "gat/{sample}.spp.nucleotide-overlap.counts",
        spp2_nuc = "gat/{sample}.spp2.nucleotide-overlap.counts",
        macs2_nuc = "gat/{sample}.macs2.nucleotide-overlap.counts",
        spp = "gat/{sample}_seg_spp",
        macs2 = "gat/{sample}_seg_macs2",
        spp2 = "gat/{sample}_seg_spp2",
    params:
        N = len(sample),
        genome = "data/dm6.chr.startstop",
        sample = "{sample}"
    wrapper:
        "file:wrappers/gat"

rule intervalstats1:
    input:
        spp = "peak_out/spp/{sample}/{sample}.bed.narrowPeak.sorted.tr",
        spp2 = "peak_out/spp2/{sample}/{sample}.bed.narrowPeak.sorted.tr",
        macs2 = "peak_out/macs2/{sample}/{sample}_peaks.narrowPeak.sorted.tr",
    output:
        "intervalstats/{sample}_query_{algo}.full",
    params:
        genome = "data/dm6.bed",
        sample = "{sample}"
    wrapper:
        "file:wrappers/intervalstats1"

rule intervalstats2:
    input:
        spp = "peak_out/spp/{sample}/{sample}.bed.narrowPeak.sorted.tr",
        spp2 = "peak_out/spp2/{sample}/{sample}.bed.narrowPeak.sorted.tr",
        macs2 = "peak_out/macs2/{sample}/{sample}_peaks.narrowPeak.sorted.tr",
        full = "intervalstats/{sample}_query_{algo}.full"
    output:
        "intervalstats/{sample}_query_{algo}.csv"
    wrapper:
        "file:wrappers/intervalstats2"

rule colocalization:
    input:
        intervalstats = expand("intervalstats/{sample}_query_{algo}.full", sample=d.keys(), algo=peakcallers),
        gat = expand("gat/{sample}.{algo}.nucleotide-overlap.counts", sample=d.keys(), algo=peakcallers)
    output:
        intervalstats = "intervalstats/intervalstast_heatmap.png",
        gat = "gat/gat_heatmap.png",
        colocalization = "colocalization/colocalization_heatmap.png"
    wrapper:
        "file:wrappers/colocalization"

rule sphinx:
    input:
        heatmap_rp=expand("heatmaps/reference_point/{algo}_{sample}_heatmap.pdf", algo=peakcallers, sample=sample),
        heatmap_sr=expand("heatmaps/scale_regions/{algo}_{sample}_heatmap.pdf",algo=peakcallers,sample=sample),
        spp_cc=expand("peak_out/spp/{sample}/{sample}.bed.pdf", sample=sample),
        macs2_model=expand("peak/macs2/{sample}/{sample}_model.pdf",sample=sample),
        fingerprint=expand("fingerprints/{sample}.pdf",sample=sample),
    output:
        "sphinx/index.rst"
    wrapper:
        "file:wrappers/sphinx"

# vim: ft=python
