rm(list = ls());tmp <- lapply(c('ggplot2','reshape2','gridExtra'), library, character.only=T)
setwd('/Users/wangz21/Library/CloudStorage/OneDrive-TheMountSinaiHospital/Huang_lab/manuscripts/ADASE/analysis/ASE_individual_analysis/ASE_Var_RegionVsClinical')

#### Load in related information
# Statistics of Sig/Sug
tmp <- as.data.frame(readxl::read_xlsx('../../../data/SampleInAnalysis/SampleInAnalysis_TenReads.xlsx'), stringsAsFactors = F)
sta <- as.data.frame(readxl::read_xlsx('../../ASE_individual/ASE_basic_SigSta_TenReads.xlsx'), stringsAsFactors = F)
sta <- sta[sta$Sampleid%in%tmp$Sampleid,]
tmp <- as.data.frame(readxl::read_xlsx('../../ROSMAP/data/SampleInAnalysis/SampleInAnalysis_TenReads.xlsx'), stringsAsFactors = F)
sta_ROSMAP <- as.data.frame(readxl::read_xlsx('../../ROSMAP/ASE_individual/ASE_basic_SigSta_TenReads.xlsx'), stringsAsFactors = F)
sta_ROSMAP <- sta_ROSMAP[sta_ROSMAP$Sampleid%in%tmp$Sampleid,]
# Sample classification
classification <- as.data.frame(readxl::read_xlsx('../../../data/SampleClassification/SampleClassification.xlsx'),stringsAsFactors=F)
classification_ROSMAP <- as.data.frame(readxl::read_xlsx('../../ROSMAP/data/SampleClassification/SampleClassification.xlsx'),stringsAsFactors=F)
# Clinical information
clinical <- read.table('../../../data/ReadCount/msbb.meta.v3.d.tsv',header = T,sep = '\t')
clinical$Classification <- classification$Classification[match(clinical$Sampleid,classification$Sampleid)]
clinical_ROSMAP <- read.table('../../ROSMAP/data/ReadCount/meta.tsv',sep = '\t',header = T,stringsAsFactors = F)
clinical_ROSMAP$Classification <- classification_ROSMAP$Classification[match(clinical_ROSMAP$Sample,classification_ROSMAP$Sample)]
colnames(clinical_ROSMAP)[match('Sample',colnames(clinical_ROSMAP))] <- 'Sampleid'

#### process data
sta_clinical <- merge(sta,clinical)
sta_clinical_ROSMAP <- merge(sta_ROSMAP,clinical_ROSMAP)

sta_clinical$SigPercent <- sta_clinical$AutoHighMAPQReadsSufSamFreqSig/sta_clinical$AutoHighMAPQReadsSufSamFreq
sta_clinical_ROSMAP$SigPercent <- sta_clinical_ROSMAP$AutoHighMAPQReadsSufSamFreqSig/sta_clinical_ROSMAP$AutoHighMAPQReadsSufSamFreq
colnames(sta_clinical_ROSMAP)[match('apoe_genotype',colnames(sta_clinical_ROSMAP))] <- 'APOE'
colnames(sta_clinical_ROSMAP)[match('pmi',colnames(sta_clinical_ROSMAP))] <- 'PMI'
colnames(sta_clinical_ROSMAP)[match('AOD',colnames(sta_clinical_ROSMAP))] <- 'AOD'
colnames(sta_clinical_ROSMAP)[match('sex',colnames(sta_clinical_ROSMAP))] <- 'SEX'
colnames(sta_clinical_ROSMAP)[match('Exonic.rate',colnames(sta_clinical_ROSMAP))] <- 'Exonic.Rate'
colnames(sta_clinical_ROSMAP)[match('rRNA.rate',colnames(sta_clinical_ROSMAP))] <- 'rRNA.rate'

ContinueousVars <- c('AOD','Exonic.Rate','PlaqueMean','PMI','rRNA.rate')
ContinueousVars_ROSMAP <- c('AOD','Exonic.Rate','PMI','rRNA.rate')
apply(sta_clinical[,ContinueousVars],2,class)
apply(sta_clinical[,ContinueousVars],2,function(x){sum(is.na(x))})
apply(sta_clinical[,ContinueousVars],2,median)
apply(sta_clinical[,ContinueousVars],2,function(x){sum(x>median(x))})
apply(sta_clinical[,ContinueousVars],2,function(x){sum(x<=median(x))})
apply(sta_clinical[,ContinueousVars],2,function(x){sum(x>=median(x))})
apply(sta_clinical[,ContinueousVars],2,function(x){sum(x<median(x))})
for (i in 1:length(ContinueousVars)) {
  tmp <- rep('Higher',nrow(sta_clinical))
  tmpValue <- sta_clinical[,ContinueousVars[i]]
  tmpMedian <- median(tmpValue)
  tmp[tmpValue<=tmpMedian] <- 'Lower'
  sta_clinical[,paste0(ContinueousVars[i],'1')] <- tmp
  
  if(ContinueousVars[i]%in%ContinueousVars_ROSMAP){
    tmp <- rep('Higher',nrow(sta_clinical_ROSMAP))
    tmpValue <- sta_clinical_ROSMAP[,ContinueousVars[i]]
    tmp[tmpValue<=tmpMedian] <- 'Lower'
    sta_clinical_ROSMAP[,paste0(ContinueousVars[i],'1')] <- tmp
  }# if
}# i

sta_clinical[is.na(sta_clinical)] <- 'NA'
sta_clinical$APOE <- factor(sta_clinical$APOE,levels = c('NA','Ambiguous','e2/e2','e3/e3','e4/e4','e2/e3','e2/e4','e3/e4'))
sta_clinical$RACE <- factor(sta_clinical$RACE,levels = c('W','H','B','A','U'))
sta_clinical$SEX <- factor(sta_clinical$SEX,levels = c('F','M'))
sta_clinical$AOD1 <- factor(sta_clinical$AOD1,levels = c('Lower','Higher'))
sta_clinical$Exonic.Rate1 <- factor(sta_clinical$Exonic.Rate1,levels = c('Lower','Higher'))
sta_clinical$PlaqueMean1 <- factor(sta_clinical$PlaqueMean1,levels = c('Lower','Higher'))
sta_clinical$PMI1 <- factor(sta_clinical$PMI1,levels = c('Lower','Higher'))
sta_clinical$rRNA.rate1 <- factor(sta_clinical$rRNA.rate1,levels = c('Lower','Higher'))

sta_clinical_ROSMAP$APOE[is.na(sta_clinical_ROSMAP$APOE)] <- 'NA'
sta_clinical_ROSMAP$APOE[sta_clinical_ROSMAP$APOE=='22'] <- 'e2/e2'
sta_clinical_ROSMAP$APOE[sta_clinical_ROSMAP$APOE=='23'] <- 'e2/e3'
sta_clinical_ROSMAP$APOE[sta_clinical_ROSMAP$APOE=='24'] <- 'e2/e4'
sta_clinical_ROSMAP$APOE[sta_clinical_ROSMAP$APOE=='33'] <- 'e3/e3'
sta_clinical_ROSMAP$APOE[sta_clinical_ROSMAP$APOE=='34'] <- 'e3/e4'
sta_clinical_ROSMAP$APOE[sta_clinical_ROSMAP$APOE=='44'] <- 'e4/e4'
sta_clinical_ROSMAP$APOE <- factor(sta_clinical_ROSMAP$APOE,levels = c('NA','Ambiguous','e2/e2','e3/e3','e4/e4','e2/e3','e2/e4','e3/e4'))
sta_clinical_ROSMAP$SEX[sta_clinical_ROSMAP$SEX=='female'] <- 'F'
sta_clinical_ROSMAP$SEX[sta_clinical_ROSMAP$SEX=='male'] <- 'M'
sta_clinical_ROSMAP$SEX <- factor(sta_clinical_ROSMAP$SEX,levels = c('F','M'))

#### boxplot
m <- rbind(c('AOD1: Higher Vs Lower','BM_10'),
           c('AOD1: Higher Vs Lower','BM_22'),
           c('AOD1: Higher Vs Lower','BM_36'),
           c('AOD1: Higher Vs Lower','DLPFC'),
           c('APOE: e2/e3 Vs e3/e3','BM_44'),
           c('RACE: H Vs B','BM_22'),
           c('RACE: W Vs B','BM_10'),
           c('RACE: W Vs B','BM_36'),
           c('RACE: W Vs H','BM_22'),
           c('SEX: F Vs M','DLPFC'))
m1 <- matrix(nrow = 0,ncol = 4); lev <- c()
for (i in 1:nrow(m)) {
  if(m[i,2]=='DLPFC'){
    tmpsta_clinical <- sta_clinical_ROSMAP
  }else{
    tmpsta_clinical <- sta_clinical[sta_clinical$Region==m[i,2],]
  }# if
  tmp1 <- unlist(strsplit(m[i,1],': '));  tmp2 <- unlist(strsplit(tmp1[2],' '))
  tmp11 <- tmp1[1]; tmp12 <- tmp2[1]; tmp13 <- tmp2[3]
  
  tmpsta_clinical_AD_Exposed <- tmpsta_clinical[tmpsta_clinical$Classification=='AD' & tmpsta_clinical[,tmp11]==tmp12,'SigPercent']
  tmpsta_clinical_AD_NotExposed <- tmpsta_clinical[tmpsta_clinical$Classification=='AD' & tmpsta_clinical[,tmp11]==tmp13,'SigPercent']
  tmpsta_clinical_Control_Exposed <- tmpsta_clinical[tmpsta_clinical$Classification=='Control' & tmpsta_clinical[,tmp11]==tmp12,'SigPercent']
  tmpsta_clinical_Control_NotExposed <- tmpsta_clinical[tmpsta_clinical$Classification=='Control' & tmpsta_clinical[,tmp11]==tmp13,'SigPercent']
  
  tmpm1 <- rbind(data.frame(Class=paste(m[i,1],m[i,2],'AD',tmp12,sep = ' '),SigPercent=tmpsta_clinical_AD_Exposed,Region=m[i,2],Classification='AD',stringsAsFactors = F),
                 data.frame(Class=paste(m[i,1],m[i,2],'AD',tmp13,sep = ' '),SigPercent=tmpsta_clinical_AD_NotExposed,Region=m[i,2],Classification='AD',stringsAsFactors = F),
                 data.frame(Class=paste(m[i,1],m[i,2],'Control',tmp12,sep = ' '),SigPercent=tmpsta_clinical_Control_Exposed,Region=m[i,2],Classification='Control',stringsAsFactors = F),
                 data.frame(Class=paste(m[i,1],m[i,2],'Control',tmp13,sep = ' '),SigPercent=tmpsta_clinical_Control_NotExposed,Region=m[i,2],Classification='Control',stringsAsFactors = F))
  m1 <- rbind(m1,tmpm1)
  
  lev <- c(lev,c(paste(m[i,1],m[i,2],'AD',tmp12,sep = ' '),paste(m[i,1],m[i,2],'AD',tmp13,sep = ' '),paste(m[i,1],m[i,2],'Control',tmp12,sep = ' '),paste(m[i,1],m[i,2],'Control',tmp13,sep = ' ')))
}# i

m1$Class <- factor(m1$Class,levels = (lev))
p <- ggplot(m1, aes(x=Class, y=SigPercent,fill = Classification)) + geom_boxplot(outlier.shape = NA)
p <- p +  scale_y_continuous(limits = quantile(m1$SigPercent, c(0.1, 0.95)))
p <- p + theme_bw() + theme(axis.text.x = element_text(size = 10,angle = 45,vjust = 1,hjust = 1),axis.text.y = element_text(size =10), 
                            axis.title.x = element_blank(), axis.title.y = element_blank(),
                            axis.ticks = element_blank(),
                            panel.grid.major = element_blank(), panel.grid.minor=element_blank(),
                            plot.title = element_text(hjust = 0.5,size = 22, face = 'bold'),
                            legend.text = element_text(size = 12,angle = 45),legend.title = element_text(size = 12),legend.position = 'none')
p <- p + geom_jitter(shape=16, position=position_jitter(0.2),color='black',size=1,alpha=0.2)
# p <- p + coord_flip()
p <- p + scale_fill_manual(values = list('AD'=rgb(228,124,124,maxColorValue = 255), 'Control'=rgb(159,184,214,maxColorValue = 255)))
ggsave(p,filename = 'BoxPlot.pdf',w=6.5,h=4,useDingbat=F)



