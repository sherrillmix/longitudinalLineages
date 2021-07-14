
source('functions.R')
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

#https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html
variantsOfConcern<-c("B.1.1.7","B.1.351","B.1.427","B.1.429","P.1")

cols<-tmpCols<-structure(c(grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(15), "#CCCCCC"),.Names= c("B.1", "B.1.1.519", "B.1.1.7 (Alpha)", "B.1.2", "B.1.243", "B.1.311", "B.1.351 (Beta)", "B.1.427 (Epsilon)", "B.1.429 (Epsilon)", "B.1.526 (Iota)", "B.1.575", "B.1.617.2 (Delta)", "B.1.621", "P.1 (Gamma)", "R.1", "Other"))
names(cols)<-sub(' .*','',names(cols))



dat<-read.csv('2021-07-07_delaware_valley_sars-cov-2_variants.csv',stringsAsFactors=FALSE)
dat$rdate<-lubridate::mdy(dat$sample_date)
#baseDate<-min(dat$rdate)
baseDate<-lubridate::mdy('4/5/2020')
dat$week<-as.numeric((dat$rdate-baseDate)/7)
dat$weekId<-1+floor(dat$week)
lineageTab<-table(dat$pango_lineage)
abundantLineages<-names(lineageTab)[lineageTab>10&names(lineageTab) %in% names(cols)]
dat$lineage<-ifelse(dat$pango_lineage %in% abundantLineages,dat$pango_lineage,'Other')
#dat$source<-ifelse(dat$rationale %in% c('asymptomatic','random'),'random',ifelse(dat$rationale=='hospitalized','hospitalized',ifelse(grepl('vaccine breakthrough',dat$rationale),'vaccine',ifelse(dat$rationale=='s drop','s drop','other'))))
dat$source<-ifelse(dat$rationale %in% c('asymptomatic','random','hospitalized','random; s drop'),'random',ifelse(grepl('vaccine breakthrough',dat$rationale),'vaccine',ifelse(dat$rationale=='s drop','s drop','other')))
countTab<-tapply(rep(1,nrow(dat)),list(dat$lineage,factor(dat$weekId,levels=1:max(dat$weekId)),dat$source),sum)
countTab[is.na(countTab)]<-0
dat$zip2<-substring(sub(' +','',sprintf('%05d',as.numeric(sub('[^0-9]+','',sub('-.*','',dat$zip))))),1,2)

mod <- rstan::stan_model("model.stan")
if(!exists('stan_sample'))stan_sample<-runStan(countTab[,,'random'],mod)
pdf('bayesTestNewData.pdf',width=6,height=6)
  plotStan(stan_sample,countTab[,,'random'],baseDate,cols)
dev.off()



mod2 <- rstan::stan_model("model_source.stan")
#runStan2<-function(countTab,hospitalTab,dropTab,vaccineTab,mod,iter=2000){
if(!exists('stan_sample2'))stan_sample2<-runStan2(countTab[,,'random'],countTab[,,'s drop'],countTab[,,'vaccine'],mod2,30000)
meanLowUp<-calcCredInt(stan_sample2,rownames(countTab))
dense<-calcDensity(stan_sample2,rownames(countTab))

pdf('bayesTestProportions.pdf',width=7,height=10)
  par(mar=c(3,3,.5,3),tcl=-.2,mgp=c(2,.3,0))
  layout(matrix(c(1,0,2,0,3,0,4,0,5,0),ncol=1),height=c(1,.01,1,.1,1,.1,1,.1,1,.2))
  plotMeanUpperLower(meanLowUp[[1]],meanLowUp[[2]],meanLowUp[[3]],countTab[,,'random'],baseDate,cols)
  par(mar=c(3,3,2,3))
  plotBars(meanLowUp$mean,baseDate,cols,showLegend=FALSE,showSum=FALSE)
  title(main='Predicted proportions data',adj=.45)
  plotBars(countTab[,,'random'],baseDate,cols,showLegend=FALSE)
  title(main='"Random" population',adj=.45,line=1)
  plotBars(countTab[,,'s drop'],baseDate,cols,showLegend=FALSE)
  title(main='S dropout',adj=.45,line=1)
  plotBars(countTab[,,'vaccine'],baseDate,cols,showLegend=FALSE)
  title(main='Vaccine breakthrough',adj=.45,line=1)
  legend('bottomleft',rownames(countTab),fill=cols[rownames(countTab)],inset=c(.1,-.6),ncol=4,xpd=NA,cex=.8)
dev.off()
pdf('bayesTestIndividualVariants.pdf',width=8,height=6)
  nCol<-4
  nRow<-ceiling(nrow(countTab)/nCol)
  layout(cbind(0,rbind(0,matrix(1:(nRow*nCol),ncol=4,byrow=TRUE),0),0),width=c(.25,rep(1,nCol),.15),height=c(.05,rep(1,nCol),.3))
  par(mar=c(0.0001,0,0,0),tcl=-.2,mgp=c(2,.3,0))
  plotIndivStan(meanLowUp[[1]],meanLowUp[[2]],meanLowUp[[3]],countTab[,,'random'],baseDate,cols)
dev.off()
pdf('bayesTestEffects.pdf',width=8,height=7)
  nCol<-4
  nRow<-ceiling(nrow(countTab)/nCol)
  layout(cbind(0,rbind(0,matrix(1:(nRow*nCol),ncol=4,byrow=TRUE),0),0),width=c(.12,rep(1,nCol),.01),height=c(.05,rep(1,nCol),.23))
  par(mar=c(2,1,1,2),tcl=-.2,mgp=c(2,.3,0))
  plotIndivDense(dense$drop,cols,xlim=c(1e-3,1e4),ylab='Fold change in odds of appearing in S1 dropout set')
  layout(cbind(0,rbind(0,matrix(1:(nRow*nCol),ncol=4,byrow=TRUE),0),0),width=c(.12,rep(1,nCol),.01),height=c(.05,rep(1,nCol),.23))
  plotIndivDense(dense$vaccine,cols,xlim=c(1e-2,1e2),ylab='Fold change in odds of appearing in vaccine breakthrough set')
dev.off()



table(substring(,1,2))
sort(unique(dat$zip))
