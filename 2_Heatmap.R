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

#### Comparison between different levels of clinical variable for AD/AsymAD/Control across brain regions
VarsMultipleGroups <- c('APOE','RACE')
VarsTwoGroups <- c('AOD1','Exonic.Rate1','PlaqueMean1','PMI1','rRNA.rate1','SEX')
m_Vars1 <- as.data.frame(matrix(nrow = length(VarsTwoGroups),ncol = 3,dimnames = list(NULL,c('variable','Exposed','NotExposed'))),stringsAsFactors=F)
m_Vars1$variable <- VarsTwoGroups; m_Vars1$Exposed <- 'Higher'; m_Vars1$NotExposed <- 'Lower'
m_Vars1$Exposed[m_Vars1$variable=='SEX'] <- 'F'; m_Vars1$NotExposed[m_Vars1$variable=='SEX'] <- 'M'
m_Vars2 <- as.data.frame(matrix(nrow = 4,ncol = 3,dimnames = list(NULL,c('variable','Exposed','NotExposed'))),stringsAsFactors=F)
m_Vars2$variable <- 'APOE'; m_Vars2$Exposed <- c('e2/e2','e2/e3','e3/e4','e4/e4'); m_Vars2$NotExposed <- 'e3/e3'
m_Vars3 <- as.data.frame(matrix(nrow = 3,ncol = 3,dimnames = list(NULL,c('variable','Exposed','NotExposed'))),stringsAsFactors=F)
m_Vars3$variable <- 'RACE'; m_Vars3$Exposed <- c('W','W','H'); m_Vars3$NotExposed <- c('H','B','B')
m_Vars <- rbind(m_Vars1,m_Vars2,m_Vars3)

Regions <- c('BM_10','BM_22','BM_36','BM_44','DLPFC'); i=1
m_fisher <- cbind(m_Vars,as.data.frame(matrix(nrow = nrow(m_Vars),ncol = length(Regions)*4,
                                              dimnames = list(NULL,paste(Regions,rep(c('Control_Pvalue','AD_Pvalue','Control_FoldChange','AD_FoldChange'),each=length(Regions)),sep = '_'))),stringsAsFactors=F))
for (i in 1:length(Regions)) {
  if(Regions[i]!='DLPFC'){
    tmpsta_clinical <- sta_clinical[sta_clinical$Region==Regions[i],]
    tmpsta_clinical_AD <- tmpsta_clinical[tmpsta_clinical$Classification=='AD',]
    tmpsta_clinical_Control <- tmpsta_clinical[tmpsta_clinical$Classification=='Control',]
  }else if(Regions[i]=='DLPFC'){
    tmpsta_clinical <- sta_clinical_ROSMAP
    tmpsta_clinical_AD <- tmpsta_clinical[tmpsta_clinical$Classification=='AD',]
    tmpsta_clinical_Control <- tmpsta_clinical[tmpsta_clinical$Classification=='Control',]
  }
  for (j in 1:nrow(m_fisher)) {
    if(m_fisher$variable[j]%in%colnames(tmpsta_clinical_AD)){
      tmpsta_clinical_Control <- tmpsta_clinical_Control[!is.na(tmpsta_clinical_Control[,m_fisher$variable[j]]),]
      tmpsta_clinical_AD <- tmpsta_clinical_AD[!is.na(tmpsta_clinical_AD[,m_fisher$variable[j]]),]
      if(sum(tmpsta_clinical_Control[,m_fisher$variable[j]]==m_fisher$Exposed[j])>=5 &
         sum(tmpsta_clinical_Control[,m_fisher$variable[j]]==m_fisher$NotExposed[j])>=5){
        m_fisher[j,paste(Regions[i],'Control_Pvalue',sep = '_')] <- wilcox.test(tmpsta_clinical_Control$SigPercent[tmpsta_clinical_Control[,m_fisher$variable[j]]==m_fisher$Exposed[j]],
                                                                                tmpsta_clinical_Control$SigPercent[tmpsta_clinical_Control[,m_fisher$variable[j]]==m_fisher$NotExposed[j]])$p.value
        m_fisher[j,paste(Regions[i],'Control_FoldChange',sep = '_')] <- median(tmpsta_clinical_Control$SigPercent[tmpsta_clinical_Control[,m_fisher$variable[j]]==m_fisher$Exposed[j]])/median(tmpsta_clinical_Control$SigPercent[tmpsta_clinical_Control[,m_fisher$variable[j]]==m_fisher$NotExposed[j]])
      }# if
      if(sum(tmpsta_clinical_AD[,m_fisher$variable[j]]==m_fisher$Exposed[j])>=5 &
         sum(tmpsta_clinical_AD[,m_fisher$variable[j]]==m_fisher$NotExposed[j])>=5){
        m_fisher[j,paste(Regions[i],'AD_Pvalue',sep = '_')] <- wilcox.test(tmpsta_clinical_AD$SigPercent[tmpsta_clinical_AD[,m_fisher$variable[j]]==m_fisher$Exposed[j]],
                                                                           tmpsta_clinical_AD$SigPercent[tmpsta_clinical_AD[,m_fisher$variable[j]]==m_fisher$NotExposed[j]])$p.value
        m_fisher[j,paste(Regions[i],'AD_FoldChange',sep = '_')] <- median(tmpsta_clinical_AD$SigPercent[tmpsta_clinical_AD[,m_fisher$variable[j]]==m_fisher$Exposed[j]])/median(tmpsta_clinical_AD$SigPercent[tmpsta_clinical_AD[,m_fisher$variable[j]]==m_fisher$NotExposed[j]])
      }# if
    }# if
  }# j
}# i
writexl::write_xlsx(m_fisher,'HeatMap.xlsx',format_headers = F)

### Plot
tmp <- paste0(m_fisher$variable,': ',m_fisher$Exposed,' Vs ',m_fisher$NotExposed)
m_fisher_plot <- expand.grid(tmp,c('BM_10_Control','BM_10_AD','BM_22_Control','BM_22_AD','BM_36_Control','BM_36_AD','BM_44_Control','BM_44_AD','DLPFC_Control','DLPFC_AD'),stringsAsFactors = F); colnames(m_fisher_plot) <- c('Clinical','Cohort')
m_fisher_plot$Pvalue <- m_fisher_plot$FoldChange <- NA
for (i in 1:nrow(m_fisher_plot)) {
  m_fisher_plot[i,'FoldChange'] <- m_fisher[match(m_fisher_plot$Clinical[i],tmp),paste0(m_fisher_plot$Cohort[i],'_FoldChange')]
  m_fisher_plot[i,'Pvalue'] <- m_fisher[match(m_fisher_plot$Clinical[i],tmp),paste0(m_fisher_plot$Cohort[i],'_Pvalue')]
}# i

m_fisher_plot$LogFC <- log2(m_fisher_plot$FoldChange)
m_fisher_plot$LogFC[is.na(m_fisher_plot$LogFC)] <- 0
m_fisher_plot$Pvalue[is.na(m_fisher_plot$Pvalue)] <- 1
m_fisher_plot$Clinical <- factor(m_fisher_plot$Clinical,levels = rev(sort(tmp)))
m_fisher_plot <- m_fisher_plot[!(m_fisher_plot$Clinical%in%c('rRNA.rate1: Higher Vs Lower','Exonic.Rate1: Higher Vs Lower','PMI1: Higher Vs Lower','PlaqueMean1: Higher Vs Lower')),]

cohorts <- unique(m_fisher_plot$Cohort);m_fisher_plot$fdr1 <- NA
for (s in 1:length(cohorts)) {ind <- which(m_fisher_plot$Cohort==cohorts[s]);m_fisher_plot$fdr1[ind] <- p.adjust(m_fisher_plot$Pvalue[ind],method = 'BH')}# for s
sum(m_fisher_plot$fdr1<0.15,na.rm = T)

p <- ggplot(data = m_fisher_plot)+geom_tile(aes(x=Cohort,y=Clinical,fill=LogFC),color=NA) + labs(title = 'ASE variant percentage comparison between different levels of variables')
p <- p + geom_tile(aes(x=Cohort,y=Clinical),fill=NA,color=ifelse(m_fisher_plot$Pvalue<0.15 & m_fisher_plot$Pvalue>=0.05,rgb(220,220,220,maxColorValue = 255),NA),size=2)
p <- p + geom_tile(aes(x=Cohort,y=Clinical),fill=NA,color=ifelse(m_fisher_plot$Pvalue<0.05 & m_fisher_plot$Pvalue>=0.01,rgb(150,150,150,maxColorValue = 255),NA),size=2)
p <- p + geom_tile(aes(x=Cohort,y=Clinical),fill=NA,color=ifelse(m_fisher_plot$Pvalue<0.01,rgb(0,0,0,maxColorValue = 255),NA),size=2)
p <- p + theme_minimal()+theme(axis.text.x = element_text(size = 25,angle = 45,vjust = 1,hjust = 1),axis.text.y = element_text(size = 25),
                               axis.title.x = element_blank(),axis.title.y = element_blank(),
                               panel.grid.major = element_blank(), plot.title = element_text(hjust = 0.5, face = 'bold',size = 25),
                               legend.text = element_text(size = 20,angle = 45,vjust = 1,hjust = 1),legend.title = element_text(size = 25),legend.position = 'bottom')
p<-p+scale_fill_gradient2(low = rgb(159,184,214,maxColorValue = 255),mid = 'white',high = rgb(228,124,124,maxColorValue = 255),midpoint = 0)
p<-p+coord_fixed(ratio = 0.5)
p<-p+geom_text(aes(x=Cohort,y=Clinical,label=ifelse(!is.na(FoldChange),round(LogFC,2),NA)),color='black',size=8)
ggsave(p,filename = 'HeatMap.pdf',w=22,h=8,useDingbat=F)


