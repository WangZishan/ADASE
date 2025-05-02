rm(list = ls())
setwd('/Users/wangz21/Library/CloudStorage/OneDrive-TheMountSinaiHospital/Huang_lab/manuscripts/ADASE')
tmp <- lapply(c('ggplot2','gridExtra','writexl','readxl','gg.layers','ggrepel','reshape2'), library, character.only=T)

#### Load in related information
### AD GWAS gene from kuan
tmp <- as.data.frame(readxl::read_xlsx('../../Huang_lab_data/GENE_LISTS/ADGWASgenes_Bellenguez_Wightman.xlsx'),stringsAsFactors=F)
Kuan_Bellenguez <- tmp$Bellenguez; Kuan_Wightman <- tmp$Wightman; Kuan_Wightman <- Kuan_Wightman[!is.na(Kuan_Wightman)]
Kuan_GWAS <- sort(unique(c(Kuan_Bellenguez,Kuan_Wightman)))
Kuan_GenesRequired <- c('APP','PSEN1','PSEN2','APOE','MAPT','SORL1','TREM2')
### AD GWAS gene set
# AD GWAS gene set Bellenguez
tmp <- as.data.frame(read_xlsx('data/ADGWASGenes/ADgwas_Bellenguez_NatGenet_2022/KnownLoci.xlsx'),stringsAsFactors=F)
GWAS_Bellenguez_Knowloci_KnownLocus <- unique(unlist(strsplit(tmp$`Known locus`,'/')))
GWAS_Bellenguez_Knowloci_ClosestGene <- unique(unlist(strsplit(tmp$Gene,'/')))
tmp <- as.data.frame(read_xlsx('data/ADGWASGenes/ADgwas_Bellenguez_NatGenet_2022/NewLoci1.xlsx'),stringsAsFactors=F)
GWAS_Bellenguez_Newloci_ClosestGene <- unique(unlist(strsplit(tmp$Gene,'/')))
tmp <- as.data.frame(read_xlsx('data/ADGWASGenes/ADgwas_Bellenguez_NatGenet_2022/NewLoci2.xlsx'),stringsAsFactors=F)
GWAS_Bellenguez_Newloci_Tier1 <- unique(tmp$Gene[tmp$`Gene Prioritization Tier`=='Tier 1'])
GWAS_Bellenguez_Newloci_Tier2 <- unique(tmp$Gene[tmp$`Gene Prioritization Tier`=='Tier 2'])

table(GWAS_Bellenguez_Knowloci_KnownLocus%in%Kuan_Bellenguez)
table(GWAS_Bellenguez_Knowloci_ClosestGene%in%Kuan_Bellenguez)
table(GWAS_Bellenguez_Newloci_ClosestGene%in%Kuan_Bellenguez)
table(GWAS_Bellenguez_Newloci_Tier1%in%Kuan_Bellenguez)
table(GWAS_Bellenguez_Newloci_Tier2%in%Kuan_Bellenguez)
# AD GWAS gene set Wightman
tmp <- as.data.frame(read_xlsx('data/ADGWASGenes/ADgwas1M_Wightman_NatGenet_2021/GWAS_Variant.xlsx'),stringsAsFactors=F)
GWAS_Wightman <- unique(unlist(strsplit(tmp$Gene,'/')))
# AD GWAS gene set other
GWAS_Bellenguez <- unique(c(GWAS_Bellenguez_Knowloci_KnownLocus,GWAS_Bellenguez_Knowloci_ClosestGene,GWAS_Bellenguez_Newloci_ClosestGene,GWAS_Bellenguez_Newloci_Tier1,GWAS_Bellenguez_Newloci_Tier2))
GWAS_Bellenguez_Responsible <- unique(c(GWAS_Bellenguez_Knowloci_KnownLocus,GWAS_Bellenguez_Newloci_Tier1))
GWAS_Bellenguez_ResponsibleIncludingTier2 <- unique(c(GWAS_Bellenguez_Knowloci_KnownLocus,GWAS_Bellenguez_Newloci_Tier1,GWAS_Bellenguez_Newloci_Tier2))
GWAS <- unique(c(GWAS_Wightman,GWAS_Bellenguez))
GWAS_Responsible <- unique(c(GWAS_Wightman,GWAS_Bellenguez_Responsible))
GWAS_ResponsibleIncludingTier2 <- unique(c(GWAS_Wightman,GWAS_Bellenguez_ResponsibleIncludingTier2))
ADGeneSet <- c('GWAS','GWAS_Responsible','GWAS_ResponsibleIncludingTier2','GWAS_Bellenguez','GWAS_Bellenguez_Responsible','GWAS_Bellenguez_ResponsibleIncludingTier2')

Kuan_GenesRequired
(Kuan_GenesRequired%in%GWAS)
(Kuan_GenesRequired%in%Kuan_GWAS)
table(Kuan_GWAS%in%GWAS); setdiff(Kuan_GWAS,GWAS)

GWAS_All <- sort(c(GWAS,Kuan_GenesRequired))
### Variant to genes
tmpfile <- gzfile('data/Statistic/VariantList2Gene/VariantList2Gene.txt.gz','r'); Var2Gene_MSBB <- read.table(tmpfile,header = T,sep = '\t',stringsAsFactors = F); close(tmpfile)
tmpfile <- gzfile('analysis/ROSMAP/data/Statistic/VariantList2Gene/VariantList2Gene.txt.gz','r'); Var2Gene_ROSMAP <- read.table(tmpfile,header = T,sep = '\t',stringsAsFactors = F); close(tmpfile)
Var2Gene <- unique(rbind(Var2Gene_MSBB,Var2Gene_ROSMAP))
### Variant HGVSc/HGVSp
tmpfile <- gzfile('data/VariantAnnotation/VariantAnnotationFinal_MSBB.txt.gz','r'); Var2HGVS_MSBB <- read.table(tmpfile,header = T,sep = '\t',stringsAsFactors = F); close(tmpfile)
tmpfile <- gzfile('data/VariantAnnotation/VariantAnnotationFinal_ROSMAP.txt.gz','r'); Var2HGVS_ROSMAP <- read.table(tmpfile,header = T,sep = '\t',stringsAsFactors = F); close(tmpfile)
Var2HGVS <- unique(rbind(Var2HGVS_MSBB,Var2HGVS_ROSMAP))
### AD Samples in analysis across brain regions
Sample_MSBB <- as.data.frame(readxl::read_xlsx('data/SampleInAnalysis/SampleInAnalysis_TenReads.xlsx'),stringsAsFactors=F)
Sample_ROSMAP <- as.data.frame(readxl::read_xlsx('analysis/ROSMAP/data/SampleInAnalysis/SampleInAnalysis_TenReads.xlsx'),stringsAsFactors=F)
Sample_ROSMAP$Region <- 'DLPFC'
Sample_All <- rbind(Sample_MSBB[,c('Region','Classification','Sampleid')],Sample_ROSMAP[,c('Region','Classification','Sampleid')])
### Load in ASE information of samples across regions
Regions <- c('BM_10','BM_22','BM_36','BM_44','DLPFC'); i=1
Statuses <- c('AD','Control'); j=1; k=1
for (i in 1:length(Regions)) {
  l <- list()
  for (j in 1:length(Statuses)) {
    tmpSample_All <- Sample_All[Sample_All$Region==Regions[i] & Sample_All$Classification==Statuses[j],]
    if(Regions[i]=='DLPFC'){tmpdir <- '../../../../OneDrive-TheMountSinaiHospital/Desktop/Subject/ADASE/ROSMAP/ASE_individual/ASE_basic_Sig'; tmpsufix <- '.ASE.WGS_het_clean.counts.tsv'
    }else{tmpdir <- '../../../../OneDrive-TheMountSinaiHospital/Desktop/Subject/ADASE/ASE_individual/ASE_basic_Sig'; tmpsufix <- '.ASE.WES_WGS_het_clean.counts.tsv'
    }# if
    
    for (k in 1:nrow(tmpSample_All)) {
      tmpASE <- read.table(file.path(tmpdir,paste0('TenReads_',tmpSample_All$Sampleid[k],tmpsufix)),header = T,sep = '\t',stringsAsFactors = F)
      if(nrow(tmpASE)>0){
        tmpASE$mapped_gene_name <- Var2Gene[match(tmpASE$HGVSg,Var2Gene$HGVSg),'mapped_gene_name']
        tmpASE1 <- tmpASE[tmpASE$mapped_gene_name%in%GWAS_All,]
        if(nrow(tmpASE1)>0){
          tmpASE1$Sampleid <- tmpSample_All$Sampleid[k];      tmpASE1$Region <- Regions[i];      tmpASE1$Status <- Statuses[j]
          l[[paste(Regions[i],Statuses[j],tmpSample_All$Sampleid[k],sep = '_')]] <- tmpASE1
        }
      }
    }# k
  }# j
  l1 <- do.call(rbind,l)
  
  l2 <- unique(l1[,c('HGVSg','Status')])
  l2 <- data.frame(l2,refCount=NA,altCount=NA,Count=NA,Count_Increase=NA,Count_Decrease=NA,Freq=NA,Freq_Increase=NA,Freq_Decrease=NA,mapped_gene_name=NA,HGVSc=NA,HGVSp=NA,Label=F,Label_HGVS=NA,stringsAsFactors = F)
  l2$mapped_gene_name <- Var2Gene[match(l2$HGVSg,Var2Gene$HGVSg),'mapped_gene_name']
  l2$HGVSc <- Var2HGVS[match(l2$HGVSg,Var2HGVS$HGVSg),'HGVSc_unique']
  l2$HGVSp <- Var2HGVS[match(l2$HGVSg,Var2HGVS$HGVSg),'HGVSp_unique']
  for (k in 1:nrow(l2)) {
    tmpl1 <- l1[l1$HGVSg==l2$HGVSg[k] & l1$Status==l2$Status[k],]
    l2$refCount[k] <- mean(tmpl1$refCount)
    l2$altCount[k] <- mean(tmpl1$altCount)
    l2$Count[k] <- nrow(tmpl1)
    l2$Count_Decrease[k] <- sum(tmpl1$altCount<tmpl1$refCount)
    l2$Count_Increase[k] <- l2$Count[k]-l2$Count_Decrease[k]
    l2$Freq[k] <- l2$Count[k]/sum(Sample_All$Region==Regions[i] & Sample_All$Classification==l2$Status[k])
    l2$Freq_Decrease[k] <- l2$Count_Decrease[k]/sum(Sample_All$Region==Regions[i] & Sample_All$Classification==l2$Status[k])
    l2$Freq_Increase[k] <- l2$Count_Increase[k]/sum(Sample_All$Region==Regions[i] & Sample_All$Classification==l2$Status[k])
    if(!is.na(l2$HGVSp[k])){
      l2$Label_HGVS[k] <- paste(l2$mapped_gene_name[k],l2$HGVSp[k],sep = ':')
    }else if((!is.na(l2$HGVSc[k])) & (is.na(l2$HGVSp[k]))){
      l2$Label_HGVS[k] <- paste(l2$mapped_gene_name[k],l2$HGVSc[k],sep = ':')
    }else if((is.na(l2$HGVSc[k])) & (is.na(l2$HGVSp[k]))){
      l2$Label_HGVS[k] <- paste(l2$mapped_gene_name[k],l2$HGVSg[k],sep = ':')
    }
  }# k
  tmpl2 <- l2[l2$Status=='AD',];  tmpl2 <- tmpl2[order(tmpl2$Freq,decreasing = T),]
  l2$Label[l2$HGVSg%in%tmpl2$HGVSg[intersect(1:5,which(tmpl2$Count>3))]] <- T
  l2$class <- 'Increase';  l2$class[l2$altCount<l2$refCount] <- 'Decrease'
  
  l3 <- l2[l2$Label,]
  write.table(l2,file.path('analysis/ASE_individual_analysis/AD_Gene_Overlap_With_Ours/GeneOverlap/',paste0(Regions[i],'_All.txt')),row.names = F,col.names = F,sep = '\t',quote = F)
  write.table(l3[l3$Label,],file.path('analysis/ASE_individual_analysis/AD_Gene_Overlap_With_Ours/GeneOverlap/',paste0(Regions[i],'_Label.txt')),row.names = F,col.names = T,sep = '\t',quote = F)
  write.table(l3$HGVSg[l3$Label],file.path('analysis/ASE_individual_analysis/AD_Gene_Overlap_With_Ours/GeneOverlap/',paste0(Regions[i],'_LabelHGVSg.txt')),row.names = F,col.names = F,sep = '\t',quote = F)
  
  Sample_Group<-c('Control'=rgb(122,183,116,maxColorValue = 255),'AD'=rgb(193,131,237,maxColorValue = 255))
  
  p <- ggplot() + geom_point(data=l2,aes(x=log10(refCount+1),y=log10(altCount+1)),alpha=0) + xlim(0,3) + ylim(0,3)
  p <- p + stat_density_2d(data=l2,aes(x=log10(refCount+1),y=log10(altCount+1),fill = after_stat(level)), geom = "polygon")
  p <- p + scale_fill_gradient2(low='white',high='black')
  p <- p + geom_abline(slope = 1,linetype = "dashed")
  p <- p + geom_point(data = l3,aes(x=log10(refCount+1),y=log10(altCount+1),color=Status,size=Freq),alpha=0.8,stroke=0) + scale_size(limits = c(0, 0.5))
  p <- p + geom_label_repel(data = l3,aes(x=log10(refCount+1),y=log10(altCount+1),label=Label_HGVS,col=Status),show.legend =FALSE,size=5)
  p <- p + scale_color_manual(name='Sample Group',values = Sample_Group)
  p <- p + labs(x='Log10 (Ref read count + 1)',y='Log10 (Alt read count + 1)',title = paste0('ASE of AD related genes in ',Regions[i],' samples'))
  p <- p + theme_bw() + theme(axis.text = element_text(size = 15),axis.title = element_text(size = 18),plot.title = element_text(hjust = 0.5, face = 'bold'),
                              legend.text = element_text(size = 15),legend.title = element_text(size = 18),legend.position = 'bottom',
                              panel.grid = element_blank())
  # p
  ggsave(p,filename=file.path('analysis/ASE_individual_analysis/AD_Gene_Overlap_With_Ours/GeneOverlap/',paste0(Regions[i],'.pdf')),h=5.5,w=5,useDingbat=F)
  
  p <- ggplot() + geom_point(data=l2,aes(x=log10(refCount+1),y=log10(altCount+1)),alpha=0) + xlim(0,3) + ylim(0,3)
  p <- p + geom_hex(data=l2,aes(x=log10(refCount+1),y=log10(altCount+1)), bins=500,binwidth = 0.15)
  p <- p + scale_fill_gradient2(low='white',high='black',limits = c(0, 20))
  p <- p + geom_abline(slope = 1,linetype = "dashed")
  p <- p + geom_point(data = l3,aes(x=log10(refCount+1),y=log10(altCount+1),color=Status,size=Freq),alpha=0.8,stroke=0) + scale_size(limits = c(0, 0.5))
  p <- p + geom_label_repel(data = l3,aes(x=log10(refCount+1),y=log10(altCount+1),label=Label_HGVS,col=Status),show.legend =FALSE,size=5)
  p <- p + scale_color_manual(name='Sample Group',values = Sample_Group)
  p <- p + labs(x='Log10 (Ref read count + 1)',y='Log10 (Alt read count + 1)',title = paste0('ASE of AD related genes in ',Regions[i],' samples'))
  p <- p + theme_bw() + theme(axis.text = element_text(size = 15),axis.title = element_text(size = 18),plot.title = element_text(hjust = 0.5, face = 'bold'),
                              legend.text = element_text(size = 15),legend.title = element_text(size = 18),legend.position = 'bottom',
                              panel.grid = element_blank())
  p
  ggsave(file.path('analysis/ASE_individual_analysis/AD_Gene_Overlap_With_Ours/GeneOverlap/',paste0(Regions[i],'_1.pdf')),h=5.5,w=5,useDingbat=F)
  
  # p <- ggplot() + geom_point(data=l2,aes(x=log10(refCount+1),y=log10(altCount+1)),alpha=0)
  # # p <- ggplot() + geom_point(data=l2,aes(x=log10(refCount+1),y=log10(altCount+1)),alpha=0.5,size=6,stroke=0,col='grey40')
  # # p <- p + stat_density_2d(data=l2,aes(x=log10(refCount+1),y=log10(altCount+1),fill = after_stat(level)), geom = "polygon")
  # p <- p + geom_hex(data=l2,aes(x=log10(refCount+1),y=log10(altCount+1)), bins=500,binwidth = 0.15)
  # p <- p + scale_fill_gradient2(low='white',high='black',limits = c(0, 50))
  # p <- p + geom_abline(slope = 1,linetype = "dashed")
  # p <- p + geom_point(data = l3,aes(x=log10(refCount+1),y=log10(altCount+1),color=Status,size=Freq),alpha=0.8,stroke=0)
  # p <- p + geom_label_repel(data = l3,aes(x=log10(refCount+1),y=log10(altCount+1),label=Label_HGVS,col=Status),show.legend =FALSE,size=4.5)
  # p <- p + scale_color_manual(name='Sample Group',values = Sample_Group)
  # p <- p + labs(x='Log10 (Ref read count + 1)',y='Log10 (Alt read count + 1)',title = paste0('ASE of AD related genes in ',Regions[i],' samples'))
  # p <- p + theme_bw() + theme(axis.text = element_text(size = 17),axis.title = element_text(size = 20),plot.title = element_text(hjust = 0.5, face = 'bold'),
  #                             legend.text = element_text(size = 17),legend.title = element_text(size = 20),
  #                             panel.grid = element_blank())
  # p
  # ggsave(file.path('analysis/ASE_individual_analysis/AD_Gene_Overlap_With_Ours/GeneOverlap/',paste0(Regions[i],'_2.pdf')),h=5,w=8,useDingbat=F)
  
  l_plot <- list()
  Genes <- sort(unique(l3$mapped_gene_name))
  for (k in 1:length(Genes)) {
    tmpl3 <- l3[l3$mapped_gene_name==Genes[k],]
    tmpl3_1 <- melt(tmpl3[,c('Label_HGVS','Status','Freq_Increase','Freq_Decrease')],id.vars = c('Label_HGVS','Status'))
    tmpl3_1$Groups <- paste(tmpl3_1$Label_HGVS,tmpl3_1$Status,sep = '*')
    
    Color<-c('Freq_Increase'=rgb(228,124,124,maxColorValue = 255),'Freq_Decrease'=rgb(159,184,214,maxColorValue = 255))
    p <- ggplot() + geom_bar(data = tmpl3_1,aes(x=Groups,y=value,fill=variable),stat = 'identity',position="stack",col='grey40') + labs(x = 'variant', y = 'Freq', title = Genes[k])
    p <- p + ylim(0,0.5)
    p <- p + theme_bw() + theme(axis.text.y = element_text(size = 20),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1,size = 20),axis.ticks.length = unit(0,'cm'),
                                axis.title =  element_blank(),plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
                                legend.title = element_text(size = 23),legend.text = element_text(size = 20),legend.position = 'none',panel.grid = element_blank())
    p <- p + scale_fill_manual(values = Color)
    # p <- p + scale_y_continuous(limits = c(0, max(unlist(l3[,c('Freq_Increase','Freq_Decrease')]))*1.1))
    p
    # ggsave(file.path('analysis/ASE_individual_analysis/AD_Gene_Overlap_With_Ours/GeneOverlap/',paste0(Regions[i],'-',Genes[k],'.pdf')),h=6,w=3,useDingbat=F)
    l_plot[[Genes[k]]] <- p
  }# k
  p1 <- grid.arrange(grobs=l_plot, nrow = 1,respect=T)
  ggsave(p1,filename = file.path('analysis/ASE_individual_analysis/AD_Gene_Overlap_With_Ours/GeneOverlap/',paste0(Regions[i],'-AllGenes.pdf')),h=10,w=20,useDingbat=F)
  
  l_plot1 <- list()
  Genes <- sort(unique(l3$mapped_gene_name))
  for (k in 1:length(Genes)) {
    tmpl3 <- l3[l3$mapped_gene_name==Genes[k],]
    tmpl3_1 <- melt(tmpl3[,c('Label_HGVS','Status','Count_Increase','Count_Decrease')],id.vars = c('Label_HGVS','Status'))
    tmpl3_1$Groups <- paste(tmpl3_1$Label_HGVS,tmpl3_1$Status,sep = '*')
    
    Color<-c('Count_Increase'=rgb(228,124,124,maxColorValue = 255),'Count_Decrease'=rgb(159,184,214,maxColorValue = 255))
    p <- ggplot() + geom_bar(data = tmpl3_1,aes(x=Groups,y=value,fill=variable),stat = 'identity',position="stack",col='grey40') + labs(x = 'variant', y = 'Freq', title = Genes[k])
    # p <- p + ylim(0,0.5)
    p <- p + theme_bw() + theme(axis.text.y = element_text(size = 20),axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1,size = 20),axis.ticks.length = unit(0,'cm'),
                                axis.title =  element_blank(),plot.title = element_text(hjust = 0.5,size = 18, face = 'bold'),
                                legend.title = element_text(size = 23),legend.text = element_text(size = 20),legend.position = 'none',panel.grid = element_blank())
    p <- p + scale_fill_manual(values = Color)
    if(Regions[i]!='DLPFC'){
      p <- p + scale_y_continuous(limits = c(0, 30))
    }else(
      p <- p + scale_y_continuous(limits = c(0, max(unlist(l3[,c('Count_Increase','Count_Decrease')]))*1.1))
    )
    p
    # ggsave(file.path('analysis/ASE_individual_analysis/AD_Gene_Overlap_With_Ours/GeneOverlap/',paste0(Regions[i],'-',Genes[k],'.pdf')),h=6,w=3,useDingbat=F)
    l_plot1[[Genes[k]]] <- p
  }# k
  p1 <- grid.arrange(grobs=l_plot1, nrow = 1,respect=T)
  ggsave(p1,filename = file.path('analysis/ASE_individual_analysis/AD_Gene_Overlap_With_Ours/GeneOverlap/',paste0(Regions[i],'-AllGenes_1.pdf')),h=10,w=20,useDingbat=F)
  
  cat(max(l3$Freq)); cat('   '); cat(max(unlist(l3[,c('Freq_Increase','Freq_Decrease')]))); cat('   '); cat(max(unlist(l3[,c('Count_Increase','Count_Decrease')]))); cat('\n')
}# i


