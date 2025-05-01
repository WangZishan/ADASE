rm(list = ls()); tmp <- lapply(c('ggplot2','reshape2','gridExtra'),library,character.only=T)
setwd('/Users/wangz21/Library/CloudStorage/OneDrive-TheMountSinaiHospital/Huang_lab/manuscripts/ADASE/analysis/ASE_individual_analysis/ASE_Var_RegionVsChrBand')

# chromosome band information
ChrBand <- read.table('../../../data/ChrBand/ChrBand.tsv',header = F,sep = '\t',stringsAsFactors = F); colnames(ChrBand) <- c('Chr','Start','End','Band','Source'); ChrBand$ChrBand <- paste(ChrBand$Chr,':',ChrBand$Band,sep = '')

classification <- c('AD','Control','AsymAD'); j=1
longInfo <- matrix(nrow = 0,ncol = 1)
# MSBB cohort
for (j in 1:length(classification)) {
  # Significance of SampleVar
  significance <- read.table(paste('ASE_SigVar_RegionVsChrBand-SampleVarSignificance_',classification[j],'.txt',sep = ''),header = T,sep = '\t',stringsAsFactors = F)
  for (k in 1:ncol(significance)) {significance[,k] <- p.adjust(significance[,k],method = 'BH')}# for k
  significance[is.na(significance)] <- 1
  tmp1 <- setdiff(ChrBand$ChrBand,rownames(significance)); tmp2 <- matrix(1,nrow = length(tmp1),ncol = ncol(significance),dimnames = list(tmp1,colnames(significance)))
  significance <- rbind(significance,tmp2); significance$ChrBand <- rownames(significance)
  significance_long <- melt(significance,id.vars = 'ChrBand'); colnames(significance_long) <- c('ChrBand','Region','Significance')
  significance_long$ChrBand <- as.character(significance_long$ChrBand); significance_long$Region <- as.character(significance_long$Region); significance_long$Significance <- as.numeric(significance_long$Significance)
  
  # Count of significant SampleVar
  tmpCount <- as.data.frame(readxl::read_xlsx(file.path('ASE_Var_RegionVsChrBand',paste(classification[j],'_SampleVar.xlsx',sep = '')),sheet = 'Sig_SampleVar_Count'),stringsAsFactors=F)
  colnames(tmpCount)[grep('ChrRegion',colnames(tmpCount))] <- 'ChrBand'; tmpCount[is.na(tmpCount)] <- 0
  tmp1 <- setdiff(ChrBand$ChrBand,tmpCount$ChrBand); tmp2 <- as.data.frame(matrix(0,nrow = length(tmp1),ncol = ncol(tmpCount),dimnames = list(NULL,colnames(tmpCount))),stringsAsFactors = F); tmp2$ChrBand <- tmp1
  tmpCount <- rbind(tmpCount,tmp2)
  tmpCount_long <- melt(tmpCount,id.vars = 'ChrBand'); colnames(tmpCount_long) <- c('ChrBand','Region','Count')
  tmpCount_long$ChrBand <- as.character(tmpCount_long$ChrBand); tmpCount_long$Region <- as.character(tmpCount_long$Region); tmpCount_long$Count <- as.numeric(tmpCount_long$Count)
  # Proportion of significant SampleVar
  tmpProportion <- as.data.frame(readxl::read_xlsx(file.path('ASE_Var_RegionVsChrBand',paste(classification[j],'_SampleVar.xlsx',sep = '')),sheet = 'Sig_SampleVar_Proportion'),stringsAsFactors=F)
  colnames(tmpProportion)[grep('ChrRegion',colnames(tmpProportion))] <- 'ChrBand'; tmpProportion[is.na(tmpProportion)] <- 0
  tmp1 <- setdiff(ChrBand$ChrBand,tmpProportion$ChrBand); tmp2 <- as.data.frame(matrix(0,nrow = length(tmp1),ncol = ncol(tmpProportion),dimnames = list(NULL,colnames(tmpProportion))),stringsAsFactors = F); tmp2$ChrBand <- tmp1
  tmpProportion <- rbind(tmpProportion,tmp2)
  tmpProportion_long <- melt(tmpProportion,id.vars = 'ChrBand'); colnames(tmpProportion_long) <- c('ChrBand','Region','Proportion')
  tmpProportion_long$ChrBand <- as.character(tmpProportion_long$ChrBand); tmpProportion_long$Region <- as.character(tmpProportion_long$Region); tmpProportion_long$Proportion <- as.numeric(tmpProportion_long$Proportion)
  # Merge information of all kinds
  tmp_long <- merge(tmpCount_long,tmpProportion_long);  tmp_long <- merge(tmp_long,significance_long); tmp_long$Proportion_Label <- round(tmp_long$Proportion,4)
  tmp <- unlist(gregexpr(':',tmp_long$ChrBand)); tmp_long$Chr <- substr(tmp_long$ChrBand,1,tmp-1); tmp_long$Band <- substr(tmp_long$ChrBand,tmp+1,1000); tmp_long$classification <- classification[j]
  
  longInfo <- rbind(longInfo,tmp_long)
}# for j
# ROSMAP cohort
for (j in 1:length(classification)) {
  # Significance of SampleVar
  significance <- read.table(file.path('../../ROSMAP/ASE_individual_analysis/ASE_Var_RegionVsChrBand',paste('ASE_SigVar_RegionVsChrBand-SampleVarSignificance_',classification[j],'.txt',sep = '')),header = T,sep = '\t',stringsAsFactors = F)
  for (k in 1:ncol(significance)) {significance[,k] <- p.adjust(significance[,k],method = 'BH')}# for k
  significance[is.na(significance)] <- 1
  tmp1 <- setdiff(ChrBand$ChrBand,rownames(significance)); tmp2 <- matrix(1,nrow = length(tmp1),ncol = ncol(significance),dimnames = list(tmp1,colnames(significance)))
  significance <- rbind(significance,tmp2); significance$ChrBand <- rownames(significance)
  significance_long <- melt(significance,id.vars = 'ChrBand'); colnames(significance_long) <- c('ChrBand','Region','Significance')
  significance_long$ChrBand <- as.character(significance_long$ChrBand); significance_long$Region <- as.character(significance_long$Region); significance_long$Significance <- as.numeric(significance_long$Significance)
  
  # Count of significant SampleVar
  tmpCount <- as.data.frame(readxl::read_xlsx(file.path('../../ROSMAP/ASE_individual_analysis/ASE_Var_RegionVsChrBand/ASE_Var_RegionVsChrBand/',paste(classification[j],'_SampleVar.xlsx',sep = '')),sheet = 'Sig_SampleVar_Count'),stringsAsFactors=F)
  colnames(tmpCount)[grep('ChrRegion',colnames(tmpCount))] <- 'ChrBand'; tmpCount[is.na(tmpCount)] <- 0
  tmp1 <- setdiff(ChrBand$ChrBand,tmpCount$ChrBand); tmp2 <- as.data.frame(matrix(0,nrow = length(tmp1),ncol = ncol(tmpCount),dimnames = list(NULL,colnames(tmpCount))),stringsAsFactors = F); tmp2$ChrBand <- tmp1
  tmpCount <- rbind(tmpCount,tmp2)
  tmpCount_long <- melt(tmpCount,id.vars = 'ChrBand'); colnames(tmpCount_long) <- c('ChrBand','Region','Count')
  tmpCount_long$ChrBand <- as.character(tmpCount_long$ChrBand); tmpCount_long$Region <- as.character(tmpCount_long$Region); tmpCount_long$Count <- as.numeric(tmpCount_long$Count)
  # Proportion of significant SampleVar
  tmpProportion <- as.data.frame(readxl::read_xlsx(file.path('../../ROSMAP/ASE_individual_analysis/ASE_Var_RegionVsChrBand/ASE_Var_RegionVsChrBand/',paste(classification[j],'_SampleVar.xlsx',sep = '')),sheet = 'Sig_SampleVar_Proportion'),stringsAsFactors=F)
  colnames(tmpProportion)[grep('ChrRegion',colnames(tmpProportion))] <- 'ChrBand'; tmpProportion[is.na(tmpProportion)] <- 0
  tmp1 <- setdiff(ChrBand$ChrBand,tmpProportion$ChrBand); tmp2 <- as.data.frame(matrix(0,nrow = length(tmp1),ncol = ncol(tmpProportion),dimnames = list(NULL,colnames(tmpProportion))),stringsAsFactors = F); tmp2$ChrBand <- tmp1
  tmpProportion <- rbind(tmpProportion,tmp2)
  tmpProportion_long <- melt(tmpProportion,id.vars = 'ChrBand'); colnames(tmpProportion_long) <- c('ChrBand','Region','Proportion')
  tmpProportion_long$ChrBand <- as.character(tmpProportion_long$ChrBand); tmpProportion_long$Region <- as.character(tmpProportion_long$Region); tmpProportion_long$Proportion <- as.numeric(tmpProportion_long$Proportion)
  # Merge information of all kinds
  tmp_long <- merge(tmpCount_long,tmpProportion_long);  tmp_long <- merge(tmp_long,significance_long); tmp_long$Proportion_Label <- round(tmp_long$Proportion,4)
  tmp <- unlist(gregexpr(':',tmp_long$ChrBand)); tmp_long$Chr <- substr(tmp_long$ChrBand,1,tmp-1); tmp_long$Band <- substr(tmp_long$ChrBand,tmp+1,1000); tmp_long$classification <- classification[j]
  
  longInfo <- rbind(longInfo,tmp_long)
}# for j
longInfo1 <- longInfo[longInfo$Significance<0.15,]
tmp <- unique(longInfo1[,c('ChrBand','Chr','Band')]); dim(tmp); table(tmp$Chr)


SigChrBand_Diff <- SigChrBand_HighFreq <- SigChrBand_Other <- c()
SigChrBand_All <- unique(longInfo1$ChrBand); i=1
Regions <- c('BM_10','BM_22','BM_36','BM_44','DLPFC'); j=1
for (i in 1:length(SigChrBand_All)) {
  for (j in 1:length(Regions)) {
    tmp <- longInfo[longInfo$ChrBand==SigChrBand_All[i] & longInfo$Region==Regions[j],]
    tmpAD <- tmp$Significance[tmp$classification=='AD']; tmpControl <- tmp$Significance[tmp$classification=='Control']; tmpAsymAD <- tmp$Significance[tmp$classification=='AsymAD']
    if((tmpAD<0.01 & tmpControl==1)|(tmpAD==1 & tmpControl<0.01)){
      SigChrBand_Diff <- c(SigChrBand_Diff,SigChrBand_All[i])
    }
    if(any(tmp$Proportion>0.2)){
      SigChrBand_HighFreq <- c(SigChrBand_HighFreq,SigChrBand_All[i])
    }
  }# j
}# i
length(SigChrBand_Diff <- ChrBand$ChrBand[ChrBand$ChrBand%in%SigChrBand_Diff])
length(SigChrBand_HighFreq <- ChrBand$ChrBand[ChrBand$ChrBand%in%SigChrBand_HighFreq])


### plot for chr band of different classes
# switch SigChrBand class and output name to switch
tmpSigChrBand <- SigChrBand_Diff; outputName <- 'Diff'; plot_h=11.5
# tmpSigChrBand <- SigChrBand_HighFreq; outputName <- 'HighFreq'; plot_h=8


longInfo2 <- longInfo[longInfo$ChrBand%in%tmpSigChrBand,]
Regions <- c('BM_10','BM_22','BM_36','BM_44','DLPFC'); i=1
l_CountPlot <- l_ProportionPlot <- list()
for (i in 1:length(Regions)) {
  tmp <- longInfo2[longInfo2$Region==Regions[i],]
  tmp$ChrBand <- factor(tmp$ChrBand,levels = rev(ChrBand$ChrBand))
  p <- ggplot(data = tmp)+geom_tile(aes(x=classification,y=ChrBand,fill=Count),color=NA) + labs(title = paste(Regions[i],sep = ' '))
  p <- p + geom_tile(aes(x=classification,y=ChrBand),fill=NA,color=ifelse(tmp$Significance<0.15 & tmp$Significance>=0.05,rgb(220,220,220,maxColorValue = 255),NA),size=0.5)
  p <- p + geom_tile(aes(x=classification,y=ChrBand),fill=NA,color=ifelse(tmp$Significance<0.05 & tmp$Significance>=0.01,rgb(150,150,150,maxColorValue = 255),NA),size=0.5)
  p <- p + geom_tile(aes(x=classification,y=ChrBand),fill=NA,color=ifelse(tmp$Significance<0.01,rgb(0,0,0,maxColorValue = 255),NA),size=0.5)
  p <- p + theme_minimal()+theme(axis.text.x = element_text(size = 25,angle = 90,vjust = 0.5,hjust = 1),axis.text.y = element_text(size = 25),
                                 axis.title.x = element_blank(),axis.title.y = element_blank(),
                                 panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold',size = 25),
                                 legend.text = element_text(size = 20,angle = 45,vjust = 1,hjust = 1),legend.title = element_text(size = 25),legend.position = 'bottom')
  
  p<-p+scale_fill_gradient2(low = 'white',high = rgb(228,124,124,maxColorValue = 255),name='Count')
  p<-p+coord_fixed(ratio = 0.5)
  p<-p+geom_text(aes(x=classification,y=ChrBand,label=Count),color='black',size=8)
  l_CountPlot[[i]] <- p
  # ggsave(p,file=file.path('ASE_Var_RegionVsChrBand',paste(Regions[i],'_SampleVar_Count.pdf',sep = '')),h=40,w=8)
  
  tmp$Proportion_Label1 <- round(tmp$Proportion*100,2)
  p <- ggplot(data = tmp)+geom_tile(aes(x=classification,y=ChrBand,fill=Proportion),color=NA) + labs(title = paste(Regions[i],sep = ' '))
  p <- p + geom_tile(aes(x=classification,y=ChrBand),fill=NA,color=ifelse(tmp$Significance<0.15 & tmp$Significance>=0.05,rgb(220,220,220,maxColorValue = 255),NA),size=0.5)
  p <- p + geom_tile(aes(x=classification,y=ChrBand),fill=NA,color=ifelse(tmp$Significance<0.05 & tmp$Significance>=0.01,rgb(150,150,150,maxColorValue = 255),NA),size=0.5)
  p <- p + geom_tile(aes(x=classification,y=ChrBand),fill=NA,color=ifelse(tmp$Significance<0.01,rgb(0,0,0,maxColorValue = 255),NA),size=0.5)
  p <- p + theme_minimal()+theme(axis.text.x = element_text(size = 25,angle = 90,vjust = 0.5,hjust = 1),axis.text.y = element_text(size = 25),
                                 axis.title.x = element_blank(),axis.title.y = element_blank(),
                                 panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold',size = 25),
                                 legend.text = element_text(size = 20,angle = 45,vjust = 1,hjust = 1),legend.title = element_text(size = 25),legend.position = 'bottom')
  if(outputName == 'Diff'){
    p<-p+scale_fill_gradient2(low = 'white',high = rgb(228,124,124,maxColorValue = 255),name='Frequency [%]',limits=c(0,0.15))
  }else{
    p<-p+scale_fill_gradient2(low = 'white',high = rgb(228,124,124,maxColorValue = 255),name='Frequency [%]',limits=c(0,1))
  }
  p<-p+coord_fixed(ratio = 0.5)
  p<-p+geom_text(aes(x=classification,y=ChrBand,label=Proportion_Label1),color='black',size=8)
  l_ProportionPlot[[i]] <- p
  # ggsave(p,file=file.path('ASE_Var_RegionVsChrBand',paste(Regions[i],'_SampleVar_Proportion.pdf',sep = '')),h=40,w=8)
}# i
p1 <- grid.arrange(grobs=l_CountPlot, nrow = 1)
ggsave(p1,filename = file.path('ASE_Var_RegionVsChrBand',paste('SampleVar_Count_',outputName,'.pdf',sep = '')),h=plot_h,w=40,useDingbat=F)
p1 <- grid.arrange(grobs=l_ProportionPlot, nrow = 1)
ggsave(p1,filename = file.path('ASE_Var_RegionVsChrBand',paste('SampleVar_Proportion_',outputName,'.pdf',sep = '')),h=plot_h,w=40,useDingbat=F)

