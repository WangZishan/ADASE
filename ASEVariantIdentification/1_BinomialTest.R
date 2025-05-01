rm(list = ls())
indir <- '/sc/arion/projects/adineto/ase'
outdir<-'/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/ADASE/ASE_individual/ASE_basic/'

### ASE analysis of each sample-variant by using binom.test
files <- setdiff(dir(indir),c('msbb.meta.tsv','msbb.meta.v2.tsv','ROSMAP'))
for (i in 1:length(files)) {
  allele_count<-read.table(file.path(indir,files[i]),header = T,sep = '\t',stringsAsFactors = F)
  p <- rep(NA,nrow(allele_count))
  for (j in 1:nrow(allele_count)) {p[j]<-binom.test(c(allele_count[j,'refCount'],allele_count[j,'altCount']),p=0.5,alternative = 'two.sided')$p.value}# for j
  p1 <- p <- unlist(p); allele_count$pval <- p
  p[!((allele_count$refCount+allele_count$altCount)>=20)] <- NA; allele_count$fdr <- p.adjust(p,method = 'BH')
  p1[!(((allele_count$refCount+allele_count$altCount)>=20)&(allele_count$refCount>=2)&(allele_count$altCount>=2))] <- NA; allele_count$fdr1 <- p.adjust(p1,method = 'BH')
  write.table(allele_count,file.path(outdir,files[i]),row.names = F,col.names = T,quote = F,sep = '\t')
}# for i

