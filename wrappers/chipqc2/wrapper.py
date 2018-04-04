from snakemake import shell
import csv

# Run chipqc on the input using the experiment-wise approach
shell('module load R; \
        Rscript /data/Lei_student/Hussain/ChipSeq/ChIPQCrds.R \
        -i=samples.csv \
        -o={snakemake.output.metadata}_spp \
        -c chr2L chr2R chr3L chr3R')

shell('module load R; \
        Rscript /data/Lei_student/Hussain/ChipSeq/ChIPQCrds.R \
        -i=samples1.csv \
        -o={snakemake.output.metadata}_macs2 \
        -c chr2L chr2R chr3L chr3R')

# Markdown generator of ChipQCAnalysis in stitch-mode
shell('module load R; \
       Rscript /data/Lei_student/Hussain/ChipSeq/ChIPQCAnalysis.R \
        -i=samples2.csv \
        -o={snakemake.output.st_html} \
        -c chr2L chr2R chr3L chr3R')

shell('module load R;  \
        Rscript /data/Lei_student/Hussain/ChipSeq/ChIPQCAnalysis.R \
        -i=samples3.csv \
        -o={snakemake.output.st_html} \
        -c chr2L chr2R chr3L chr3R')

# Markdown geneartor of ChipQCAnalysis in experiment-mode
shell('module load R; \
        Rscript /data/Lei_student/Hussain/ChipSeq/ChIPQCAnalysis.R \
        -i={snakemake.input.metadata} \
        -o={snakemake.output.exp_html} \
        -c chr2L chr2R chr3L chr3R')
