from snakemake import shell
import pandas

_df = pandas.read_table(str(snakemake.input.full), names=['query', 'closest_ref', 'length', 'distance', 'numerator', 'denominator', 'pval'])
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

df = pandas.DataFrame([dict(
            query=snakemake.input.spp,
            filename=str(snakemake.input.full),
            reference=snakemake.input.macs2,
            n=float(n),
            n_05=n_05,
            n_01=n_01,
            n_001=n_001,
            f_05=f_05,
            f_01=f_01,
            f_001=f_001)]
)
df.tocsv(str("intervalstats/{snakemake.output}", sep='\t', index=False)

