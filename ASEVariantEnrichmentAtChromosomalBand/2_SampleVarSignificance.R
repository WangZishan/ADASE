rm(list = ls()); setwd('/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/ADASE/')

# All sample in analysis clinical information
clinicals <- as.data.frame(readxl::read_xlsx('data/SampleInAnalysis/SampleInAnalysis_TenReads.xlsx'),stringsAsFactors=F)
#clinicals <- clinicals[clinicals$Sampleid%in%c('B18C014.hB_RNA_8825','B18C014.hB_RNA_8835','B18C014.hB_RNA_8845','B18C014.hB_RNA_8865'),]
# variant to chromosome band
infile <- gzfile('data/Statistic/VariantList2Gene/VariantList2Gene_MappedToChrBand.txt.gz','r'); Var2ChrBand <- read.table(infile,header = T,sep = '\t',stringsAsFactors = F); close(infile)
Var2ChrBand$ChrBand <- paste('chr',Var2ChrBand$Chromsome,':',Var2ChrBand$mapped_chr_region,sep = '')

# Permutation analysis
classification <- c('AD','AsymAD','Control'); regions <- c('BM_10','BM_22','BM_36','BM_44')
for (k in 1:length(classification)) {
  clinical <- clinicals[clinicals$Classification==classification[k],]
  Count <- as.data.frame(readxl::read_xlsx(file.path('ASE_individual_analysis/ASE_Var_RegionVsChrBand/ASE_Var_RegionVsChrBand/',paste(classification[k],'_SampleVar.xlsx',sep = '')),sheet = 'Sig_SampleVar_Count'),stringsAsFactors=F); rownames(Count) <- Count[,'ChrRegion']; Count <- Count[,setdiff(colnames(Count),'ChrRegion')]
  
  m <- matrix(0,nrow = nrow(Count),ncol = ncol(Count),dimnames = list(rownames(Count),colnames(Count)))
  for (j in 1:length(regions)) {
    tmpclinical <- clinical[clinical$Region==regions[j],]
    all <- list()
    for (i in 1:nrow(tmpclinical)) {
      sig <- read.table(file.path('ASE_individual/ASE_basic_Sig',paste('TenReads_',tmpclinical$Sampleid[i],'.ASE.WES_WGS_het_clean.counts.tsv',sep = '')),header = T,sep = '\t',stringsAsFactors = F)
      notsig <- read.table(file.path('ASE_individual/ASE_basic_NotSig',paste('TenReads_',tmpclinical$Sampleid[i],'.ASE.WES_WGS_het_clean.counts.tsv',sep = '')),header = T,sep = '\t',stringsAsFactors = F)
      
      sig$class <- 'Significant'; notsig$class <- 'Not'; tmpall <- rbind(sig,notsig)
      tmpall$HGVSg <- paste(substr(tmpall$contig,4,nchar(tmpall$contig)),':g.',tmpall$position,tmpall$refAllele,'>',tmpall$altAllele,sep = ''); tmpall$ChrBand <- Var2ChrBand$ChrBand[match(tmpall$HGVSg,Var2ChrBand$HGVSg)]
      all[[i]] <- tmpall[,c('ChrBand','class')]
    }# for i
    all <- do.call(rbind,all)
    
    rand<-100000; tmpCount <- Count[,regions[j]]
    for (i in 1:rand) {
      tmpall <- all; tmpall$class <- sample(all$class,nrow(all))
      tmpall <- tmpall[tmpall$class=='Significant',]
      tmpall_count<- table(tmpall$ChrBand)[rownames(Count)]; tmpall_count[is.na(tmpall_count)] <- 0
      m[,regions[j]] <- m[,regions[j]]+as.numeric(tmpall_count>=tmpCount)
    }# for i
  }# for j
  m <- m/rand
  write.table(m,file.path('ASE_individual_analysis/ASE_Var_RegionVsChrBand',paste('ASE_SigVar_RegionVsChrBand-SampleVarSignificance_',classification[k],'.txt',sep = '')),row.names = T,col.names = T,quote = F,sep = '\t')
}# for k

