rm(list = ls()); setwd('/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/ADASE/')
outdir <- 'ASE_individual_analysis/ASE_Var_RegionVsChrBand/ASE_Var_RegionVsChrBand'

# All sample in analysis clinical information
clinicals <- as.data.frame(readxl::read_xlsx('data/SampleInAnalysis/SampleInAnalysis_TenReads.xlsx'),stringsAsFactors=F)
#clinicals <- clinicals[clinicals$Sampleid%in%c('B18C014.hB_RNA_8825','B18C014.hB_RNA_8835','B18C014.hB_RNA_8845','B18C014.hB_RNA_8865'),]
# variant to chromosome band
infile <- gzfile('data/Statistic/VariantList2Gene/VariantList2Gene_MappedToChrBand.txt.gz','r'); Var2ChrBand <- read.table(infile,header = T,sep = '\t',stringsAsFactors = F); close(infile)
Var2ChrBand$ChrBand <- paste('chr',Var2ChrBand$Chromsome,':',Var2ChrBand$mapped_chr_region,sep = '')
# chromosome band
ChrBand <- read.table('data/ChrBand/ChrBand.tsv',header = F,sep = '\t',stringsAsFactors = F); ChrBand$ChrBand <- paste(ChrBand$V1,':',ChrBand$V4,sep = ''); ChrBand <- ChrBand$ChrBand[ChrBand$ChrBand%in%Var2ChrBand$ChrBand]

# SampleVariant among chromosome band vs regions
regions <- c('BM_10','BM_22','BM_36','BM_44'); classification <- c('AD','AsymAD','Control')
for (r in 1:length(classification)) {
  clinical <- clinicals[clinicals$Classification==classification[r],]
  
  # require variables
  Vars <- c('SampleVar','Sample','Var','SamplePerVar','VarPerSample')
  for (i in 1:length(Vars)) {
    tmp <- matrix(nrow = length(ChrBand),ncol = length(regions),dimnames = list(ChrBand,regions))
    assign(paste('m_',Vars[i],'_All',sep = ''),tmp); assign(paste('m_',Vars[i],'_Sig',sep = ''),tmp)
  }# for
  # Count
  for (i in 1:length(regions)) {
    tmpclinical <- clinical[clinical$Region==regions[i],]
    tmplist_all <- tmplist_sig <- list()
    for (j in 1:nrow(tmpclinical)) {
      tmpsig <- read.table(file.path('ASE_individual/ASE_basic_Sig/',paste('TenReads_',tmpclinical$Sampleid[j],'.ASE.WES_WGS_het_clean.counts.tsv',sep = '')),sep = '\t',stringsAsFactors = F,header = T)
      tmpnot <- read.table(file.path('ASE_individual/ASE_basic_NotSig/',paste('TenReads_',tmpclinical$Sampleid[j],'.ASE.WES_WGS_het_clean.counts.tsv',sep = '')),sep = '\t',stringsAsFactors = F,header = T)
	  if(class(tmpsig$refAllele)=='logical'){tmpsig$refAllele <- as.character(rep('T',nrow(tmpsig)))}
      if(class(tmpsig$altAllele)=='logical'){tmpsig$altAllele <- as.character(rep('T',nrow(tmpsig)))}
      if(class(tmpnot$refAllele)=='logical'){tmpnot$refAllele <- as.character(rep('T',nrow(tmpnot)))}
      if(class(tmpnot$altAllele)=='logical'){tmpnot$altAllele <- as.character(rep('T',nrow(tmpnot)))}
      tmp <- rbind(tmpsig,tmpnot)
      tmplist_sig[[j]] <- data.frame(Sample=tmpclinical$Sampleid[j],tmpsig,stringsAsFactors = F)
      tmplist_all[[j]] <- data.frame(Sample=tmpclinical$Sampleid[j],tmp,stringsAsFactors = F)
    }# for
    tmplist_sig <- do.call(rbind,tmplist_sig); tmplist_sig$HGVSg <- paste(substr(tmplist_sig$contig,4,nchar(tmplist_sig$contig)),':g.',tmplist_sig$position,tmplist_sig$refAllele,'>',tmplist_sig$altAllele,sep = ''); tmplist_sig$ChrBand <- Var2ChrBand$ChrBand[match(tmplist_sig$HGVSg,Var2ChrBand$HGVSg)]
    tmplist_all <- do.call(rbind,tmplist_all); tmplist_all$HGVSg <- paste(substr(tmplist_all$contig,4,nchar(tmplist_all$contig)),':g.',tmplist_all$position,tmplist_all$refAllele,'>',tmplist_all$altAllele,sep = ''); tmplist_all$ChrBand <- Var2ChrBand$ChrBand[match(tmplist_all$HGVSg,Var2ChrBand$HGVSg)]
    
    for (k in 1:length(ChrBand)) {m_SampleVar_Sig[ChrBand[k],regions[i]] <- sum(tmplist_sig$ChrBand==ChrBand[k])}# for
    for (k in 1:length(ChrBand)) {m_Sample_Sig[ChrBand[k],regions[i]] <- length(unique(tmplist_sig$Sample[tmplist_sig$ChrBand==ChrBand[k]]))}# for
    for (k in 1:length(ChrBand)) {m_Var_Sig[ChrBand[k],regions[i]] <- length(unique(tmplist_sig$HGVSg[tmplist_sig$ChrBand==ChrBand[k]]))}# for
    
    for (k in 1:length(ChrBand)) {m_SampleVar_All[ChrBand[k],regions[i]] <- sum(tmplist_all$ChrBand==ChrBand[k])}# for
    for (k in 1:length(ChrBand)) {m_Sample_All[ChrBand[k],regions[i]] <- length(unique(tmplist_all$Sample[tmplist_all$ChrBand==ChrBand[k]]))}# for
    for (k in 1:length(ChrBand)) {m_Var_All[ChrBand[k],regions[i]] <- length(unique(tmplist_all$HGVSg[tmplist_all$ChrBand==ChrBand[k]]))}# for
  }# for i
  m_SamplePerVar_All <- m_SampleVar_All/m_Var_All; m_VarPerSample_All <- m_SampleVar_All/m_Sample_All
  m_SamplePerVar_Sig <- m_SampleVar_Sig/m_Var_Sig; m_VarPerSample_Sig <- m_SampleVar_Sig/m_Sample_Sig
  
  # Output the above count and Significant Variant Proportion
  for (i in 1:length(Vars)) {
    tmp1 <- get(paste('m_',Vars[i],'_All',sep = '')); tmp2 <- get(paste('m_',Vars[i],'_Sig',sep = '')); tmp3 <- tmp2/tmp1
    l <- list()
    l[[paste('All_',Vars[i],'_Count',sep = '')]] <- data.frame(ChrRegion=rownames(tmp1),tmp1,stringsAsFactors = F)
    l[[paste('Sig_',Vars[i],'_Count',sep = '')]] <- data.frame(ChrRegion=rownames(tmp2),tmp2,stringsAsFactors = F)
    l[[paste('Sig_',Vars[i],'_Proportion',sep = '')]] <- data.frame(ChrRegion=rownames(tmp3),tmp3,stringsAsFactors = F)
    writexl::write_xlsx(l,file.path(outdir,paste(classification[r],'_',Vars[i],'.xlsx',sep = '')),format_headers = T)
  }# for i
}# for r

