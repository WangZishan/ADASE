rm(list = ls())
tmp <- lapply(c('ggplot2','reshape2'),library,character.only=T)
indir <- 'C:/Users/wangz21/OneDrive - The Mount Sinai Hospital/Desktop/Subject/ADASE'
setwd('C:/Users/wangz21/Box Sync/Huang_lab/manuscripts/ADASE')

### Load related information
# variant mapped to gene for variants
tmpfile <- gzfile('data/Statistic/VariantList2Gene/VariantList2ProteinCodingGene.txt.gz','r') ;var2gene1 <- read.table(tmpfile,header = T,sep = '\t',stringsAsFactors = F); close(tmpfile)
tmpfile <- gzfile('ROSMAP/data/Statistic/VariantList2Gene/VariantList2ProteinCodingGene.txt.gz','r') ;var2gene2 <- read.table(tmpfile,header = T,sep = '\t',stringsAsFactors = F); close(tmpfile)
var2gene <- unique(rbind(var2gene1,var2gene2))
# required matrix
Regions <- c('BM_10','BM_22','BM_36','BM_44')
m <- matrix(nrow = length(Regions),ncol = 10,dimnames = list(Regions,c('DLPFC_All','DLPFC_SampleCutoff','DLPFC_FDR0.15','MSBB_All','MSBB_SampleCutoff','MSBB_Pval0.05','Intersect_All','Intersect_PositiveConsistent','Intersect_NegativeConsistent','Intersect_NonConsistent')))
# DLPFC ASE information
ASE <- read.table(file.path(indir,'ROSMAP/ASE_group/ASE_basic/ASE_basic_DLPFC.txt'),header = T,sep = '\t',stringsAsFactors = F); m[,'DLPFC_All'] <- nrow(ASE)
ASE <- ASE[(ASE$ADSampleCount > 5) & (ASE$ControlSampleCount > 5) & (ASE$ADSamplePercent > 0) & (ASE$ControlSamplePercent > 0),]; m[,'DLPFC_SampleCutoff'] <- nrow(ASE)
ASE$FDR <- p.adjust(ASE$InteractionPvalue,method = 'BH')
DLPFC <- ASE[(ASE$FDR < 0.15) & (!is.na(ASE$FDR)),]; m[,'DLPFC_FDR0.15'] <- nrow(DLPFC)
### significant variant: MSBB Pval < 0.05 + DLPFC FDR < 0.15
for (j in 1:length(Regions)) {
  ASE <- read.table(file.path(indir,'ASE_group/ASE_basic',paste('ASE_basic_',Regions[j],'.txt',sep = '')),header = T,sep = '\t',stringsAsFactors = F); m[j,'MSBB_All'] <- nrow(ASE)
  ASE <- ASE[(ASE$ADSampleCount > 5) & (ASE$ControlSampleCount > 5) & (ASE$ADSamplePercent > 0) & (ASE$ControlSamplePercent > 0),]; m[j,'MSBB_SampleCutoff'] <- nrow(ASE)
  ASE$FDR <- p.adjust(ASE$InteractionPvalue,method = 'BH')
  MSBB <- ASE[(ASE$InteractionPvalue < 0.05) & (!is.na(ASE$InteractionPvalue)),]; m[j,'MSBB_Pval0.05'] <- nrow(MSBB)
  
  Inter_All <- intersect(MSBB$HGVSg,DLPFC$HGVSg); m[j,'Intersect_All'] <- length(Inter_All)
  Inter_Positive <- intersect(MSBB$HGVSg[MSBB$InteractionEstimate > 0],DLPFC$HGVSg[DLPFC$InteractionEstimate > 0]); m[j,'Intersect_PositiveConsistent'] <- length(Inter_Positive)
  Inter_Negative <- intersect(MSBB$HGVSg[MSBB$InteractionEstimate < 0],DLPFC$HGVSg[DLPFC$InteractionEstimate < 0]); m[j,'Intersect_NegativeConsistent'] <- length(Inter_Negative)
  Inter_Non <- setdiff(Inter_All,c(Inter_Positive,Inter_Negative)); m[j,'Intersect_NonConsistent'] <- length(Inter_Non)
  
  tmpMSBB <- MSBB; colnames(tmpMSBB) <- paste('MSBB',colnames(MSBB),sep = '_')
  tmpDLPFC <- DLPFC; colnames(tmpDLPFC) <- paste('DLPFC',colnames(DLPFC),sep = '_')
  tmp1 <- tmp2 <- tmp3 <- c()
  if(length(Inter_Positive)>0){tmp1 <- cbind(Consistency='PositiveConsistent',tmpMSBB[match(Inter_Positive,tmpMSBB$MSBB_HGVSg),],tmpDLPFC[match(Inter_Positive,tmpDLPFC$DLPFC_HGVSg),])}
  if(length(Inter_Negative)>0){tmp2 <- cbind(Consistency='NegativeConsistent',tmpMSBB[match(Inter_Negative,tmpMSBB$MSBB_HGVSg),],tmpDLPFC[match(Inter_Negative,tmpDLPFC$DLPFC_HGVSg),])}
  if(length(Inter_Non)>0){tmp3 <- cbind(Consistency='NonConsistent',tmpMSBB[match(Inter_Non,tmpMSBB$MSBB_HGVSg),],tmpDLPFC[match(Inter_Non,tmpDLPFC$DLPFC_HGVSg),])}
  output <- rbind(tmp1,tmp2,tmp3)
  output <- cbind(output,var2gene[match(output$MSBB_HGVSg,var2gene$HGVSg),setdiff(colnames(var2gene),'HGVSg')])
  write.table(output,file.path('ASE_group/ASE_basic_Sig',paste('ASE_basic_Pval0.05_',Regions[j],'_FDR0.15_DLPFC.txt',sep = '')),row.names = F,col.names = T,quote = F,sep = '\t')
}# for
writexl::write_xlsx(data.frame(Region=rownames(m),m,stringsAsFactors = F),'ASE_group/ASE_basic_Sig_3.xlsx',format_headers = F)

