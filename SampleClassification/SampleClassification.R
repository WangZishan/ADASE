rm(list = ls())
setwd('C:/Users/wangz21/Box Sync/Huang_lab/manuscripts/ADASE/')

# Load in clinical information
clinical <- read.table('data/ReadCount_DNAHetVariant/msbb.meta.d.tsv',header = T,sep = '\t',stringsAsFactors = F)

# Classification of clinical samples
clinical$Classification <- NA
for (i in 1:nrow(clinical)) {
  if(!is.na(clinical$CERAD[i]) & !is.na(clinical$bbscore[i]) & !is.na(clinical$CDR[i])){
    if(((clinical$CERAD[i]%in%c('NL','possAD') & clinical$bbscore[i]%in%c(0,1,2)) | (clinical$CERAD[i]%in%c('NL') & clinical$bbscore[i]%in%c(3))) & clinical$CDR[i]<1){
      clinical$Classification[i] <- 'Control'
    }else if(clinical$CERAD[i]%in%c('possAD','probAD','defAD') & clinical$bbscore[i]%in%c(3,4,5,6) & clinical$CDR[i]<1){
      clinical$Classification[i] <- 'AsymAD'
    }else if(clinical$CERAD[i]%in%c('probAD','defAD') & clinical$bbscore[i]%in%c(3,4,5,6) & clinical$CDR[i]>=1){
      clinical$Classification[i] <- 'AD'
    }else{
      clinical$Classification[i] <- 'Exclude'
    }
  }# for if
}# for i
clinical$Classification[is.na(clinical$Classification)] <- 'NA'
writexl::write_xlsx(clinical,'data/SampleClassification/SampleClassification.xlsx',format_headers = F)



clinical[clinical$Classification=='control' & clinical$bbscore==3 & !is.na(clinical$Classification) & !is.na(clinical$bbscore),]
table(clinical$Classification)
table(clinical$Classification,clinical$Region)

