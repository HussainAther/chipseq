from snakemake import shell
import pandas


shell(  '''
        IntervalStats \
        --query {snakemake.input.spp} \
        --reference {snakemake.input.macs2} \
        --output intervalstats/{snakemake.params.sample}_query_spp.full \
        --domain {snakemake.params.genome} \
        ''')

shell("touch {snakemake.output}")

shell(  '''
        IntervalStats \
        --query {snakemake.input.macs2} \
        --reference {snakemake.input.spp} \
        --output intervalstats/{snakemake.params.sample}_query_macs2.full \
        --domain {snakemake.params.genome} \
        ''')

shell("touch {snakemake.output}")

shell(  '''
        IntervalStats \
        --query {snakemake.input.spp2} \
        --reference {snakemake.input.macs2} \
        --output intervalstats/{snakemake.params.sample}_query_spp2.full \
        --domain {snakemake.params.genome} \
        ''')

shell("touch {snakemake.output}")
