source('functions.R')
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

#https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html
variantsOfConcern<-c("B.1.1.7","B.1.351","B.1.427","B.1.429","P.1")

cols<-structure(c("yellow2", "palegreen", "lightseagreen", "skyblue", "steelblue1", "royalblue", "navy", "purple", "hotpink", "indianred", "red", "orange","lightgoldenrod1", "honeydew4"),.Names=c("B.1", "B.1.1.7 (UK)", "B.1.2", "B.1.243", "B.1.311", "B.1.351 (South Africa)", "B.1.427 (California)", "B.1.429 (California)", "B.1.525", "B.1.526 (NY)", "B.1.617 (India)", "P.1 (Brazil)", "P.2", "Other"))
names(cols)<-sub(' .*','',names(cols))
cols<-apply(col2rgb(cols),2,function(xx)do.call(rgb,c(as.list(xx),maxColorValue=255)))


dat<-read.csv('adm_2021-05-21_confidence-intervals_delaware-valley-sars-cov-2-variants.csv',stringsAsFactors=FALSE)
dat$lineage2<-ifelse(dat$lineage=='B.1','Other',dat$lineage)
dat$rdate<-lubridate::mdy(dat$date_by_week)
baseDate<-min(dat$rdate)
dat$week<-1+(dat$rdate-baseDate)/7
dat$lineage3<-ifelse(apply(do.call(rbind,lapply(variantsOfConcern,grepl,dat$lineage2,fixed=TRUE)),2,any),dat$lineage2,'Other')

countTab<-tapply(dat$counts,list(dat$lineage3,factor(dat$week,levels=1:max(dat$week))),sum)
countTab[is.na(countTab)]<-0
countTab<-countTab[order(rownames(countTab)!='Other'),]
rownames(countTab)<-sub('Other','Not VoC',sub(' .*','',rownames(countTab)))

mod <- rstan::stan_model("model.stan")
mod2 <- rstan::stan_model("model_fix1.stan")

if(!exists('stan_sample'))stan_sample<-runStan(countTab,mod)

pdf('bayesTest.pdf',width=6,height=6)
  tmp<-cols
  names(tmp)[names(tmp)=='Other']<-'Not VoC'
  plotStan(stan_sample,countTab,baseDate,tmp)
dev.off()


countTab2<-tapply(dat$counts,list(dat$lineage2,factor(dat$week,levels=1:max(dat$week))),sum)
countTab2[is.na(countTab2)]<-0
countTab2<-countTab2[order(rownames(countTab2)!='Other'),]
countTab2<-countTab2[apply(countTab2,1,sum)>0,]
rownames(countTab2)<-sub('Other','Not VoI/VoC',sub(' .*','',rownames(countTab2)))

if(!exists('stan_sample2'))stan_sample2<-runStan(countTab2,mod)

pdf('bayesTest2.pdf',width=6,height=6)
  tmp<-cols
  names(tmp)[names(tmp)=='Other']<-'Not VoI/VoC'
  plotStan(stan_sample2,countTab2,baseDate,tmp)
dev.off()

countTab3<-tapply(dat$counts,list(dat$lineage,factor(dat$week,levels=1:max(dat$week))),sum)
countTab3[is.na(countTab3)]<-0
countTab3<-countTab3[order(rownames(countTab3)!='Other'),]
countTab3<-countTab3[apply(countTab3,1,sum)>0,]
rownames(countTab3)<-sub(' .*','',rownames(countTab3))
if(!exists('stan_sample3'))stan_sample3<-runStan(countTab3,mod)
if(!exists('stan_sample4'))stan_sample4<-runStan(countTab3,mod2,iter=10000)

pdf('bayesTest3.pdf',width=6,height=6)
  meanLowUp<-plotStan(stan_sample4,countTab3,baseDate,cols)
dev.off()
meanLowUp[[1]]
