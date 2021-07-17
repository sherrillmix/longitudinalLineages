
source('functions.R')

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
dat$zip2<-substring(sub(' +','',sprintf('%05d',as.numeric(sub('[^0-9]+','',sub('-.*','',dat$zip))))),1,2)
countTab<-tapply(rep(1,nrow(dat)),list(dat$lineage,factor(dat$weekId,levels=1:max(dat$weekId)),dat$source),sum)
countTab[is.na(countTab)]<-0

dnar::withAs(xx=dat[dat$zip2 %in% c('08','17','19'),],table(xx$lineage,xx$weekId,xx$zip2))
