if(!dir.exists('out'))dir.create('out')
source('functions.R')
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)
set.seed(12345)


#https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html
voi<-list(
  'Alpha'=c('B.1.1.7',unique(dat$Pango_Lineage[grep('^Q\\.',dat$Pango_Lineage)])),
  'Beta'=c('B.1.351',unique(dat$Pango_Lineage[grep('^B.1.351\\.',dat$Pango_Lineage)])),
  'Gamma'=c('P.1',unique(dat$Pango_Lineage[grep('^P\\.1\\.',dat$Pango_Lineage)])),
  'Epsilon'=c('B.1.427','B.1.429'),
  'Eta'='B.1.525',
  'Iota'='B.1.526',
  'Kappa'='B.1.617.1',
  'Mu'=c('B.1.621', 'B.1.621.1'),
  'Zeta'='P.2'
) 
voc<-list('Delta'=c('B.1.617.2', unique(dat$Pango_Lineage[grep('^AY\\.',dat$Pango_Lineage)])))
greekLookup<-structure(c(rep(names(voi),sapply(voi,length)),rep(names(voc),sapply(voc,length))),.Names=c(unlist(voi),unlist(voc)))

baseDate<-lubridate::mdy('3/29/2020')
dat<-read.csv('tableS1.csv',stringsAsFactors=FALSE)
dat$rdate<-lubridate::mdy(dat$Collection_Date)
dat$week<-as.numeric((dat$rdate-baseDate)/7)
dat$weekId<-1+floor(dat$week)
dat$greek<-ifelse(dat$Pango_Lineage %in% names(greekLookup),greekLookup[dat$Pango_Lineage],'Other')
dat$source<-ifelse(dat$Rationale %in% c('asymptomatic','surveillance','hospitalized'),'surveillance',dat$Rationale)
if(any(dateSelect<-is.na(dat$weekId)))warning(sum(dateSelect),' samples missing date')

lineageTab<-table(dat$Pango_Lineage)
abundantLineages<-names(lineageTab)[lineageTab>10&names(lineageTab)!='None']
dat$abundant<-ifelse(dat$Pango_Lineage %in% abundantLineages,dat$Pango_Lineage,'Other')
dat[dat$abundant=='Other'&dat$greek!='Other','abundant']<-sprintf('Other %s',dat[dat$abundant=='Other'&dat$greek!='Other','greek'])
countTab<-tapply(rep(1,nrow(dat)),list(factor(dat$abundant,levels=c('Other',sort(unique(dat$abundant[dat$abundant!='Other'])))),factor(dat$weekId,levels=1:max(dat$weekId)),dat$source),sum)
countTab[is.na(countTab)]<-0
greek<-structure(ifelse(rownames(countTab) %in% names(greekLookup),greekLookup[rownames(countTab)],'Other'),.Names=rownames(countTab))

mod <- rstan::stan_model("model_sourceVariants.stan")
if(!exists('stan_sample'))stan_sample<-runStan(countTab[,,'surveillance'],countTab[,,'s drop'],countTab[,,'vaccine breakthrough'],greek,mod,3000)

meanLowUp<-calcCredInt(stan_sample$stan,names=rownames(countTab))
dense<-calcDensity(stan_sample$stan,names=rownames(countTab),greekNames=names(stan_sample$greek))
meanLowUpVac<-calcCredInt(stan_sample$stan,names=rownames(countTab),colPattern='^propsVaccine\\[')
meanLowUpDrop<-calcCredInt(stan_sample$stan,names=rownames(countTab),colPattern='^propsDrop\\[')
dropTab<-getTableDrop(stan_sample)
strainTab<-getTable(stan_sample,FALSE)

greeks<-structure(ifelse(grepl('Other',rownames(countTab)),sub('Other ','',rownames(countTab)),ifelse(rownames(countTab) %in% names(greekLookup),greekLookup[rownames(countTab)],'Other')),.Names=rownames(countTab))
greekOrder<-rownames(countTab)[order(greeks=='Other',greeks,rownames(countTab)=='Other',rownames(countTab))]

reds<-c('#c44b07', '#fb5394', '#d17448', '#f70028', '#ed819e', '#d2003f', '#ff9669', '#95435d', '#9b5449') #https://medialab.github.io/iwanthue/
bluePurp<-c("#2090ff", "#ec2de4", "#019cdb", "#b70e8e", "#344ef8", "#a85d88", "#9fa2e0", "#ff9cef", "#72549d")
bluePurp2<-c("#db76b1", "#038bff", "#ef40b1", "#c1baff", "#d645f8", "#cb9dd8", "#7769ff", "#924883", "#684ad0", "#8e70a8", "#ab1da8", "#4757cb", "#d498ff")
cols<-c(
  B.1.1.7 = "#33a02c", `Other Alpha` = '#44b13d',
  `Other Beta` = "#a6cee3", 
  AY.12 = reds[8], AY.14 = reds[2], AY.20 = reds[3], AY.24 = reds[1], AY.25 = reds[5], AY.3 = reds[9], AY.4 = reds[7], B.1.617.2 =  reds[4], `Other Delta` = reds[6],
  B.1.525 ="#ff7f00",
  `Other Gamma` =  "#EEE685",
  P.1 = "#1CFFCE", P.1.2 = "#3EFFF1",
  B.1.526 = "#FFD700",
  `Other Kappa` = "#b2df8a",
  B.1 = bluePurp2[2], B.1.1 = bluePurp2[7], B.1.1.434 = bluePurp2[3], B.1.1.519 = bluePurp2[4], B.1.2 = bluePurp2[12], B.1.234 = bluePurp2[1], B.1.243 = bluePurp2[6], B.1.311 = bluePurp2[8], B.1.575 = bluePurp2[9], B.1.596 = bluePurp2[10], B.1.621 = bluePurp2[11], B.1.637 = bluePurp2[5], R.1 = bluePurp2[13],
  Other = "#CCCCEE"
)


pdf('out/countsPredictions.pdf',width=7,height=8)
  par(mar=c(1,2.9,.5,2.6),tcl=-.2,mgp=c(2,.3,0))
  layout(matrix(c(0,1,0,2,0,3,0,4,0),ncol=1),height=c(.15,1,.15,1,.15,1,.15,1,.62))
  plotBars(countTab[greekOrder,,'random'],baseDate,cols,showLegend=FALSE,showDates=FALSE)
  title(main='Surveillance',adj=.45,line=1,xpd=NA)
  plotBars(countTab[greekOrder,,'s drop'],baseDate,cols,showLegend=FALSE,showDates=FALSE)
  title(main='S dropout',adj=.45,line=1,xpd=NA)
  plotBars(countTab[greekOrder,,'vaccine breakthrough'],baseDate,cols,showLegend=FALSE,showDates=FALSE)
  title(main='Vaccine breakthrough',adj=.45,line=1,xpd=NA)
  plotBars(meanLowUp$mean[greekOrder,],baseDate,cols,showLegend=FALSE,showSum=FALSE,skipDate=2,dateCex=.9)
  title(main='Estimated surveillance proportions',adj=.45,line=1,xpd=NA)
  legend('bottomleft',names(cols),pt.bg=cols,inset=c(-.05,-.8),ncol=8,xpd=NA,cex=.8,pt.cex=1.8,pch=22)
dev.off()
pdf('out/sDropEnrichment.pdf',width=8,height=8)
  nCol<-5
  nRow<-ceiling(nrow(countTab)/nCol)
  layout(cbind(0,rbind(0,matrix(1:(nRow*nCol),ncol=nCol,byrow=TRUE),0),0),width=c(.12,rep(1,nCol),.01),height=c(.05,rep(1,nRow),.23))
  par(mar=c(2,1,1,2),tcl=-.2,mgp=c(2,.3,0))
  plotIndivDense(dense$drop[greekOrder],cols,xlim=c(1e-3,1e3),ylab='Fold change in odds of appearing in S1 dropout set')
dev.off()
pdf('out/vaccineEnrichment.pdf',width=8,height=8)
  nCol<-5
  nRow<-ceiling(nrow(countTab)/nCol)
  layout(cbind(0,rbind(0,matrix(1:(nRow*nCol),ncol=nCol,byrow=TRUE),0),0),width=c(.12,rep(1,nCol),.01),height=c(.05,rep(1,nRow),.23))
  par(mar=c(2,1,1,2),tcl=-.2,mgp=c(2,.3,0))
  plotIndivDense(dense$vaccine[greekOrder],cols,xlim=c(1e-2,1e2),ylab='Fold change in odds of appearing in vaccine breakthrough set')
dev.off()
pdf('out/individualVariants.pdf',width=8,height=8)
  nCol<-5
  nRow<-ceiling(nrow(countTab)/nCol)
  layout(cbind(0,rbind(0,matrix(1:(nRow*nCol),ncol=nCol,byrow=TRUE),0),0),width=c(.3,rep(1,nCol),.15),height=c(.05,rep(1,nRow),.4))
  par(mar=c(0.0001,0,0,0),tcl=-.2,mgp=c(2,.3,0))
  plotIndivStan(meanLowUp[[1]][greekOrder,],meanLowUp[[2]][greekOrder,],meanLowUp[[3]][greekOrder,],countTab[greekOrder,,'random'],baseDate,cols,nCol=nCol)
dev.off()

out<-formatC(signif(strainTab[,1:3],2))
out<-sub('^([0-9])$','\\1.0',out)
colnames(out)<-c("Mean","Lower 95% CrI","Upper 95% CrI")
out<-out[c(rownames(out)[!rownames(out) %in% greekOrder],greekOrder),]
write.csv(out,'out/strainFoldVaccineEnrichment.csv',quote=FALSE)
out<-formatC(signif(dropTab[,1:3],2))
out<-sub('^([0-9])$','\\1.0',out)
colnames(out)<-c("Mean","Lower 95% CrI","Upper 95% CrI")
write.csv(out[greekOrder,],'out/dropFoldVaccineEnrichment.csv',quote=FALSE)
out<-meanLowUp$mean
colnames(out)<-as.character(baseDate+(1:ncol(out)-1)*7)
write.csv(t(out),'out/variantProps.csv')

