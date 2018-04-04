import tempfile
from snakemake import shell

script = """
library(spp);
library(caTools)
chip.data <- read.bam.tags("{snakemake.input.ip}");
input.data <- read.bam.tags("{snakemake.input.inp}");
binding.characteristics <- get.binding.characteristics(chip.data,srange=c(50,500),bin=5,accept.all.tags=T);
print(paste("binding peak separation distance =",binding.characteristics$peak$x));
pdf(file="{snakemake.output.bed}.pdf",width=5,height=5)
par(mar=c(3.5,3.5, 1.0, 0.5), mgp=c(2,0.65,0), cex=0.8);
plot(binding.characteristics$cross.correlation,type='l',xlab="strand shift",ylab="cross-correlation");
abline(v=binding.characteristics$peak$x,lty=2,col=2)
dev.off();
chip.data.info<-select.informative.tags(chip.data,binding.characteristics);
input.data.info<-select.informative.tags(input.data,binding.characteristics);
chip.data.qua <- remove.local.tag.anomalies(chip.data.info);
input.data.qua <- remove.local.tag.anomalies(input.data.info);
tag.shift <- round(binding.characteristics$peak$x/2)
smoothed.density <- get.smoothed.tag.density(chip.data.qua,control.tags=input.data.qua,bandwidth=200,step=100,tag.shift=tag.shift);
writewig(smoothed.density,"{snakemake.output.bed}.density.wig","Example smoothed, background-subtracted tag density");
rm(smoothed.density);
smoothed.enrichment.estimate <- get.smoothed.enrichment.mle(chip.data.qua,input.data.qua,bandwidth=200,step=100,tag.shift=tag.shift)
writewig(smoothed.enrichment.estimate,"{snakemake.output.bed}.enrichment.wig","Example smoothed maximum likelihood log2 enrichment estimate");
enrichment.estimates <- get.conservative.fold.enrichment.profile(chip.data.qua,input.data.qua,fws=500,step=100,alpha=0.01);
writewig(enrichment.estimates,"{snakemake.output.narrowPeak}.enrichment.estimates.wig","Example conservative fold-enrichment/depletion estimates shown on log2 scale");
rm(enrichment.estimates);
broad.clusters <- get.broad.enrichment.clusters(chip.data.qua,input.data.qua,window.size=1e3,z.thr=3,tag.shift=round(binding.characteristics$peak$x/2))
write.broadpeak.info(broad.clusters,"{snakemake.output.bed}.broadPeak")
fdr <- 1e-2;
write.table(cbind(rep("X", length(broad.clusters$X$s)), broad.clusters$X$s, broad.clusters$X$e), file="{snakemake.output.narrowPeak}", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
detection.window.halfsize <- binding.characteristics$whs;
tec.filter = FALSE;
bp<-find.binding.positions(signal.data=chip.data.qua, control.data=input.data.qua, fdr=fdr, method=tag.lwcc, whs=detection.window.halfsize, tec.filter = FALSE)
bp <- add.broad.peak.regions(chip.data.qua,input.data.qua,bp,window.size=1000,z.thr=3)
print(paste("detected",sum(unlist(lapply(bp$npl,function(d) length(d$x)))),"peaks"));
write.narrowpeak.binding(bp,"{snakemake.output.narrowPeak}")
""".format(**locals())

script_filename = snakemake.output.bed + '.R'
fout=open(script_filename,"w")
fout.write(script)
fout.close()
shell('Rscript {script_filename} &> {snakemake.output.bed}.log')
#shell('Rscript {snakemake.output.bed}.R')
shell("cd {snakemake.params.prefix} && ln -sf $(basename {snakemake.output.narrowPeak}) $(basename {snakemake.output.bed})")
shell('touch -h {snakemake.output.narrowPeak}')
