rm(list = ls()); setwd('/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/ADASE/ASE_individual')

### Load in related information
# clinical information
clinicals <- as.data.frame(readxl::read_xlsx('../data/SampleClassification/SampleClassification.xlsx'),stringsAsFactors=F)
colnames(clinicals)[match('Sample',colnames(clinicals))] <- 'Sampleid'
clinicals <- clinicals[clinicals$Classification%in%c('AD','AsymAD','Control'),]
#clinicals <- clinicals[clinicals$Sampleid%in%c('B18C014.hB_RNA_8825','B18C014.hB_RNA_8835','B18C014.hB_RNA_8845','B18C014.hB_RNA_8865'),]
# Variant list sample frequency
infile <- gzfile('/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/ADASE/ASE_individual_analysis/ReferenceBias/ReferenceFraction/ReferenceFraction_VariantListSta_10_0.01.txt.gz','r')
VarSta <- read.table(infile,header = T,sep = '\t',stringsAsFactors = F); close(infile)
tmp <- unlist(gregexpr(':',VarSta$HGVSg)); VarSta$Chr <- as.character(substr(VarSta$HGVSg,1,tmp-1)); VarSta <- VarSta[VarSta$Chr%in%as.character(1:22),]
VarSta <- VarSta$HGVSg[VarSta$All_Freq>0.2]

m <- as.data.frame(matrix(nrow = nrow(clinicals),ncol = 8, dimnames = list(NULL,c('Sampleid','All','Auto','AutoHighMAPQ','AutoHighMAPQReadsSuf','AutoHighMAPQReadsSufSamFreq','AutoHighMAPQReadsSufSamFreqSig','AutoHighMAPQReadsSufSamFreqNotSig'))),stringsAsFactors = F)
for (i in 1:nrow(clinicals)) {
  tmpall <- read.table(file.path('ASE_basic/',paste(clinicals$Sampleid[i],'.ASE.WES_WGS_het_clean.counts.tsv',sep = '')),header = T,sep = '\t',stringsAsFactors = F)
  tmpall <- tmpall[,setdiff(colnames(tmpall),c('fdr','fdr1'))]
  tmpall$HGVSg <- paste(substr(tmpall$contig,4,nchar(tmpall$contig)),':g.',tmpall$position,tmpall$refAllele,'>',tmpall$altAllele,sep = '')
  m[i,'Sampleid'] <- clinicals$Sampleid[i];  m[i,'All'] <- nrow(tmpall)
  
  # remove variant located at XY chromosome
  tmpall <- tmpall[!(tmpall$contig%in%c('chrX','chrY')),]
  m[i,'Auto'] <- nrow(tmpall)
  
  # remove variant with low maping quality
  tmpall <- tmpall[(tmpall$lowMAPQDepth/tmpall$rawDepth)<=0.01,]
  m[i,'AutoHighMAPQ'] <- nrow(tmpall)
  
  # remove variant without enough read count
  tmpall <- tmpall[(tmpall$refCount+tmpall$altCount)>=10,]
  m[i,'AutoHighMAPQReadsSuf'] <- nrow(tmpall)
  
  # remove variant with lower frequency
  tmpall <- tmpall[tmpall$HGVSg%in%VarSta,]
  m[i,'AutoHighMAPQReadsSufSamFreq'] <- nrow(tmpall)
  
  # significant varaint
  tmpall$fdr <- p.adjust(tmpall$pval,method = 'BH')
  tmpsig <- tmpall[tmpall$fdr<0.05,]
  tmpnotsig <- tmpall[tmpall$fdr>=0.05,]
  m[i,'AutoHighMAPQReadsSufSamFreqSig'] <- nrow(tmpsig)
  m[i,'AutoHighMAPQReadsSufSamFreqNotSig'] <- nrow(tmpnotsig)
  
  write.table(tmpsig,file.path('ASE_basic_Sig',paste('TenReads_',paste(clinicals$Sampleid[i],'.ASE.WES_WGS_het_clean.counts.tsv',sep = ''),sep = '')),row.names = F,col.names = T,quote = F,sep = '\t')
  write.table(tmpnotsig,file.path('ASE_basic_NotSig',paste('TenReads_',paste(clinicals$Sampleid[i],'.ASE.WES_WGS_het_clean.counts.tsv',sep = ''),sep = '')),row.names = F,col.names = T,quote = F,sep = '\t')
}# for i
writexl::write_xlsx(m,'ASE_basic_SigSta_TenReads.xlsx',format_headers = F)
