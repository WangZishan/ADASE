rm(list = ls())

### Functions to use
glmCal = function(expr){
  InteractionEstimate <- InteractionPvalue <- AltEstimate <- AltPvalue <- CaseEstimate <- CasePvalue <- InterceptEstimate <- InterceptPvalue <- NA
  fit_res = try(withWarnings(expr), FALSE)
  if(sum(grepl("Error", fit_res))==0) {
    fit = try(fit_res$value)
    warns = sum(unlist(lapply(fit_res$warnings, function(x) {
      y=unlist(x)
      grepl("NaNs produced", y) |
        grepl("https://urldefense.proofpoint.com/v2/url?u=http-3A__glm.fit&d=DwIGAg&c=shNJtf5dKgNcPZ6Yh64b-A&r=liGlo3GLmViTeFNGvfRSNnMrcydqqtU4-S0M4YzGctE&m=yx2H_VnjxhDMsTg-LTBX8TNtbk0lf9pzWo-wCSXtAP0&s=tpHlo9bfb7pqcC66Mq3tt9x9y1y5glnmkElloXXkXBQ&e= : algorithm did not converge", y) |
        grepl("Model failed to converge",y)|
        grepl("Model is nearly unidentifiable",y)})))
    if(warns==0){
      fit_summary_coef <- summary(fit)$coefficients
      InterceptEstimate <- fit_summary_coef['(Intercept)','Estimate']; InterceptPvalue <- fit_summary_coef['(Intercept)','Pr(>|z|)']
      AltEstimate <- fit_summary_coef['Allelealt','Estimate']; AltPvalue <- fit_summary_coef['Allelealt','Pr(>|z|)']
      CaseEstimate <- fit_summary_coef['GroupAD','Estimate']; CasePvalue <- fit_summary_coef['GroupAD','Pr(>|z|)']
      InteractionEstimate <- fit_summary_coef['Allelealt:GroupAD','Estimate']; InteractionPvalue <- fit_summary_coef['Allelealt:GroupAD','Pr(>|z|)']
    }
  }
  list(InteractionEstimate=InteractionEstimate, InteractionPvalue=InteractionPvalue, AltEstimate=AltEstimate, AltPvalue=AltPvalue,
       CaseEstimate=CaseEstimate, CasePvalue=CasePvalue, InterceptEstimate=InterceptEstimate, InterceptPvalue=InterceptPvalue)
}
withWarnings <- function(expr) {
  myWarnings <- NULL
  wHandler <- function(w) {
    myWarnings <<- c(myWarnings, list(w))
    invokeRestart("muffleWarning")
  }
  val <- withCallingHandlers(expr, warning = wHandler)
  list(value = val, warnings = myWarnings)
} 

### Set up working environment
tmp <- lapply(c('readxl','lme4'),library,character.only=T)
indir <- '/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/ADASE'
setwd('/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/ADASE')

### Load in related information
# allele count information
refCount <- read.table(file.path(indir,'data/ProfileTransform/FurtherTransform/refCount_TenReads.txt'),header = T,sep = '\t',stringsAsFactors = F)
altCount <- read.table(file.path(indir,'data/ProfileTransform/FurtherTransform/altCount_TenReads.txt'),header = T,sep = '\t',stringsAsFactors = F)
# clinical information
clinical <- as.data.frame(read_xlsx('data/SampleInAnalysis/SampleInAnalysis_TenReads.xlsx'),stringsAsFactors=F)

### linear mix effect model analysis
regions <- c('BM_10','BM_22','BM_36','BM_44')
for (i in 1:length(regions)) {
  tmpsamples <- clinical$Sampleid[clinical$Region==regions[i]]
  tmprefCount <- refCount[match(tmpsamples,rownames(refCount)),]
  tmpaltCount <- altCount[match(tmpsamples,rownames(altCount)),]
  tmpclinical <- clinical[match(tmpsamples,clinical$Sampleid),]
  
  write.table(matrix(c('HGVSg','InteractionEstimate','InteractionPvalue','AltEstimate','AltPvalue','CaseEstimate','CasePvalue','InterceptEstimate','InterceptPvalue','ADSampleCount','ADSamplePercent','ControlSampleCount','ControlSamplePercent'),nrow = 1),
              file.path('ASE_group/ASE_basic/',paste('ASE_basic_',regions[i],'.txt',sep = '')),row.names = F,col.names = F,quote = F,sep = '\t')
  HGVSg <- colnames(tmprefCount); HGVSg <- substr(HGVSg,2,1000); HGVSg <- gsub('\\.g\\.',':g.',HGVSg); HGVSg <- gsub('\\.A','>A',HGVSg); HGVSg <- gsub('\\.C','>C',HGVSg); HGVSg <- gsub('\\.G','>G',HGVSg); HGVSg <- gsub('\\.T','>T',HGVSg)
  tmpADAll <- sum(tmpclinical$Classification=='AD'); tmpControlAll <- sum(tmpclinical$Classification=='Control')
  for (j in 1:ncol(tmprefCount)) {
    input <- data.frame(Count=c(tmprefCount[,j],tmpaltCount[,j]),
                        Allele=c(rep('ref',nrow(tmprefCount)),rep('alt',nrow(tmprefCount))),
                        Group=rep(tmpclinical$Classification,2),
                        Individual=rep(tmpclinical$Sampleid,2))
    input <- input[!is.na(input$Count),]
    input <- input[input$Group%in%c('AD','Control'),]
    input$Allele <- factor(input$Allele,levels = c('ref','alt'))
    input$Group <- factor(input$Group,levels = c('Control','AD'))
    
    tmpAD <- sum(input$Group=='AD')/2; tmpControl <- sum(input$Group=='Control')/2
    if((tmpAD >= 5) & (tmpControl >= 5)){
      fit = glmCal(glmer.nb(Count ~ (1|Individual) + Allele + Group + Group:Allele, data=input))
      write.table(matrix(c(HGVSg[j],fit$InteractionEstimate,fit$InteractionPvalue,fit$AltEstimate,fit$AltPvalue,fit$CaseEstimate,fit$CasePvalue,fit$InterceptEstimate,fit$InterceptPvalue,
                           tmpAD,tmpAD/tmpADAll,tmpControl,tmpControl/tmpControlAll),nrow = 1),
                  file.path('ASE_group/ASE_basic/',paste('ASE_basic_',regions[i],'.txt',sep = '')),row.names = F,col.names = F,quote = F,sep = '\t',append = T)
    }# if
  }# for j
}# for i
