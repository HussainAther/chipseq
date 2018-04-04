from snakemake import shell
import matplotlib.pyplot as plt
import pybedtools
import os

kind = snakemake.params.kind

d = {}
for bed in snakemake.input.beds:
    peakcaller = bed.split('/')[1]
    d[peakcaller] = pybedtools.BedTool(bed)

if 'found-in-all' in kind:
    anchor_peakcaller = kind.split('_')[0]
    anchor = d[anchor_peakcaller]
    for pc, bed in d.items():
        if pc == anchor:
            continue
        anchor += bed
    anchor.saveas(str(snakemake.output))
#    if anchor == "": 
#        sys.exit()
#    else:
#        anchor.saveas(str(snakemake.output))

elif '_only' in kind:
    anchor_peakcaller = kind.split('_')[0]
    anchor = d[anchor_peakcaller]
    for pc, bed in d.items():
        if pc == anchor_peakcaller:
            continue
        anchor -= bed
    anchor.saveas(str(snakemake.output))
#    if anchor == "":
#        sys.exit()
#    else:
#        anchor.saveas(str(snakemake.output))

elif '-not-' in kind:
    anchor_peakcaller, other_peakcaller = kind.split('-not-')
    anchor = d[anchor_peakcaller]
    other = d[other_peakcaller]
    result = anchor.intersect(other, v=True)
    result.saveas(str(snakemake.output))
#    if result == "":
#        sys.exit()
#    else:
#        result.saveas(str(snakemake.output))

else:
    raise ValueError('Unhandled kind: "{0}"'.format(kind))
