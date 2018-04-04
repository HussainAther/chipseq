import tempfile
from snakemake import shell

script = """
.. snakemake documentation master file, created by
sphinx-quickstart on Thu Jul 27 17:54:40 2017.
You can adapt this file completely to your liking, but it should at least
contain the root `toctree` directive.

Welcome to snakemake's documentation!
=====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:
* Peak-Caller Output
* * Reference point heatmap
* * Scale regions heatmiap
                          *          * Cross-correlational profile spp
                          *             * Model macs2
                          *                * Fingerprint
                          *
                          *
                          *
                          *                Indices and tables
                          *                ==================
                          *
                          *                * :ref:`genindex`
                          *                * :ref:`modindex`
                          *                * :ref:`search`

""".format(**locals())

script_filename = snakemake.output
fout=open(script_filename,"w")
fout.write(script)
fout.close()
with script_filename as out:
    write(script)
shell("make {script_filename} &> {snakemake.output}')
