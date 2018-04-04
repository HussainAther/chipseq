import os
from snakemake.utils import makedirs
import pandas
import yaml
shell.prefix = 'set -e; set -o pipefail; '

# The config file is expected to have the following format:
#
#   beds:
#       label1: filename1
#       label2: filename2
#       ...
#
#   domains:
#       domain_label_1: domain_filename_1
#       domain_label_2: domain_filename_2
#       ...
#
#   output: output-directory
#
#
#
config = yaml.load(open('example_config.yaml'))

# Number of shufflings for GAT
N = 1000

def from_label(wc, val):
    """
    Helper function to return the BED filename given the label
    """
    key = getattr(wc, val)
    if val == 'domain':
        return config['domains'][key]
    else:
        return config['beds'][key]


outdir = config['output']
makedirs(outdir)

targets = expand(
    '{outdir}/{algorithm}/{domain}/{query}_vs_{reference}.txt',
    outdir=outdir,
    domain=config['domains'].keys(),
    query=config['beds'].keys(),
    reference=config['beds'].keys(),
    algorithm=['IntervalStats', 'GAT', 'jaccard', 'fisher'],
)

rule target:
    input: targets


rule jaccard:
    input:
        domain=lambda wc: config['domains'][getattr(wc, 'domain')],
        query=lambda wc: config['beds'][getattr(wc, 'query')],
        reference=lambda wc: config['beds'][getattr(wc, 'reference')],
    output: '{outdir}/jaccard/{domain}/{query}_vs_{reference}.txt'
    shell:
        """
        bedtools intersect -a {input.query} -b {input.domain} | bedtools sort -i stdin > {output}.query.jaccard
        bedtools intersect -a {input.reference} -b {input.domain} | bedtools sort -i stdin > {output}.reference.jaccard
        bedtools jaccard -a {output}.query.jaccard -b {output}.reference.jaccard > {output}
        rm {output}.query.jaccard {output}.reference.jaccard
        """


rule fisher:
    input:
        domain=lambda wc: config['domains'][getattr(wc, 'domain')],
        query=lambda wc: config['beds'][getattr(wc, 'query')],
        reference=lambda wc: config['beds'][getattr(wc, 'reference')],
    output: '{outdir}/fisher/{domain}/{query}_vs_{reference}.txt'
    shell:
        """
        bedtools intersect -a {input.query} -b {input.domain} | bedtools sort -i stdin > {output}.query.fisher
        bedtools intersect -a {input.reference} -b {input.domain} | bedtools sort -i stdin > {output}.reference.fisher
        bedtools fisher -a {output}.query.fisher -b {output}.reference.fisher -g dm6.chromsizes > {output}
        rm {output}.query.fisher {output}.reference.fisher
        """


rule intervalstats:
    input:
        domain=lambda wc: config['domains'][getattr(wc, 'domain')],
        query=lambda wc: config['beds'][getattr(wc, 'query')],
        reference=lambda wc: config['beds'][getattr(wc, 'reference')],
    output: '{outdir}/IntervalStats/{domain}/{query}_vs_{reference}.txt'
    run:
        if input.query == input.reference:
            run_self = '--self'
        else:
            run_self = ''
        shell(
            '''
            IntervalStats \\
            --query {input.query} \\
            --reference {input.reference} \\
            --output {output}.full \\
            --domain {input.domain} \\
            {run_self}
            ''')

        # Summarize the output into a faster-to-parse file used by downstream
        # analysis code.
        #
        # Output has columns:
        #
        # - n_{05,01,001}: number of significant associations at {0.05, 0.01,
        #   0.001} respectively
        #
        # - f_{05,01,001}: fraction of total that are signficant
        #
        # - n: number of features
        #
        # - query, reference: labels
        #
        # - filename: "all" filename containing the details in case anything
        #   needs re-calculation.
        _df = pandas.read_table(
            str(output[0]) + '.full',
            names=['query', 'closest_ref', 'length', 'distance',
                   'numerator', 'denominator', 'pval'])

        n = float(len(_df))

        def frac(x):
            if n == 0:
                return np.nan
            return x / n

        n_05 = sum(_df.pval < 0.05)
        n_01 = sum(_df.pval < 0.01)
        n_001 = sum(_df.pval < 0.001)
        f_05 = frac(n_05)
        f_01 = frac(n_01)
        f_001 = frac(n_001)

        df = pandas.DataFrame(
            [
                dict(
                    query=wildcards.query,
                    filename=str(output[0]) + '.full',
                    reference=wildcards.reference,
                    n=float(n),
                    n_05=n_05,
                    n_01=n_01,
                    n_001=n_001,
                    f_05=f_05,
                    f_01=f_01,
                    f_001=f_001,
                )
            ]
        )
        df.to_csv(str(output[0]), sep='\t', index=False)




rule gat:
    input:
        domain=lambda wc: config['domains'][getattr(wc, 'domain')],
        query=lambda wc: config['beds'][getattr(wc, 'query')],
        reference=lambda wc: config['beds'][getattr(wc, 'reference')],
    output: '{outdir}/GAT/{domain}/{query}_vs_{reference}.txt'
    run:
        shell(
            '''
            cut -f1,2,3 {input.query} > {output}.query.tmp
            cut -f1,2,3 {input.reference} > {output}.reference.tmp
            '''
        )
        shell(
            '''
            gat-run.py \\
            --ignore-segment-tracks \\
            --annotations {output}.reference.tmp \\
            --segments {output}.query.tmp \\
            --workspace {input.domain} \\
            --counter nucleotide-overlap \\
            --cache {output}.cache \\
            --num-samples {N} \\
            --output-counts-pattern {output}.%s.counts \\
            --log {output}.log \\
            --stdout {output}
            ''')

        shell("rm {output}.query.tmp {output}.reference.tmp")



# vim: ft=python
