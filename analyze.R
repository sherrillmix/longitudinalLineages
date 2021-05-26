options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

dat<-read.csv('adm_2021-05-21_confidence-intervals_delaware-valley-sars-cov-2-variants.csv')
dat$lineage2<-ifelse(dat$lineage=='B.1','Other',dat$lineage)
dat$rdate<-lubridate::mdy(dat$date_by_week)
baseDate<-min(dat$rdate)
dat$week<-1+(dat$rdate-baseDate)/7
#https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html
variantsOfConcern<-c("B.1.1.7","B.1.351","B.1.427","B.1.429","P.1")
dat$lineage3<-ifelse(apply(do.call(rbind,lapply(variantsOfConcern,grepl,dat$lineage2,fixed=TRUE)),2,any),dat$lineage2,'Other')
countTab<-tapply(dat$counts,list(dat$lineage3,factor(dat$week,levels=1:max(dat$week))),sum)
countTab[is.na(countTab)]<-0
countTab<-countTab[order(rownames(countTab)!='Other'),]
propTab<-apply(countTab,2,function(xx)xx/sum(xx))
mod <- rstan::stan_model("model.stan")

stan_sample <- rstan::sampling(
  mod,
  data=list(
    counts=countTab,
    nTime=ncol(countTab),
    nLineage=nrow(countTab),
    firstGuess=c(4,rep(0,nrow(countTab)-1))
  ),
  iter=2000,
  chains=50,
  control=list(max_treedepth=15),
)
print(stan_sample,c('means','changes'),include=FALSE)
#print(stan_sample,c('means','changes'),include=TRUE)
pdf('test.pdf',width=30,height=30)
pairs(stan_sample,pars=c('changeSigma','means[1,1]','means[1,10]','means[2,10]')) #
dev.off()
mat<-as.matrix(stan_sample)
props<-mat[,grepl('^props',colnames(mat))]
meanCrI<-function(xx)c(mean(xx),quantile(xx,c(.025,.975)))
meanLowerUpper<-apply(props,2,meanCrI)
convertToMat<-function(xx){
  row<-as.numeric(sub('[^0-9]+([0-9]+),[0-9]+.*','\\1',names(xx)))
  col<-as.numeric(sub('[^0-9]+[0-9]+,([0-9]+).*','\\1',names(xx)))
  tapply(xx,list(row,col),c)
}
means<-convertToMat(meanLowerUpper[1,])
lower<-convertToMat(meanLowerUpper[2,])
upper<-convertToMat(meanLowerUpper[3,])

#https://sashamaps.net/docs/resources/20-colors/
cols<-c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b') 
pdf('bayesTest.pdf',width=6,height=6)
  par(mar=c(3,3,.5,2.1),tcl=-.2,mgp=c(2,.3,0))
  layout(matrix(c(1,0,2,0),ncol=1),height=c(1,.01,1,.15))
  plot(1,1,type='n',xlim=range(dat$week)+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Estimated proportion',xaxs='i')
  for(ii in 1:nrow(means)){
    lines(1:ncol(means),means[ii,],col=cols[ii])
    polygon(c(1:ncol(means),ncol(means):1),c(lower[ii,],rev(upper[ii,])),border=NA,col=sprintf('%s66',cols[ii]))
  }
  abline(h=0)
  text(ncol(means)+.5,means[,ncol(means)],sub('Other','Not VoC',sub(' .*','',rownames(countTab))),cex=.7,xpd=NA,adj=0)
  prettyDates<-lubridate::mdy(c('4/1/20','7/1/20','10/1/20','1/1/21','4/1/21'))
  axis(1,(prettyDates-baseDate)/7+1,sub('  +',' ',format(prettyDates,'%b %e %Y')),cex=.8,tcl=-.5,mgp=c(3,.7,0))
  #
  plot(1,1,type='n',xlim=range(dat$week)+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Proportion',xaxs='i')
  for(ii in 1:ncol(propTab))rect(ii-.5,cumsum(propTab[,ii]),ii+.5,cumsum(c(0,propTab[-nrow(propTab),ii])),col=cols[1:nrow(propTab)])
  uniqDate<-unique(dat[,c('week','rdate')])
  dnar::slantAxis(1,uniqDate$week,sub('  +',' ',format(uniqDate$rdate,'%b %e %Y')),location=.5,cex=.6)
  counts<-apply(countTab,2,sum)
  dnar::slantAxis(3,which(counts>0),counts[counts>0],location=.4,cex=.6,axisArgs=list(lwd=NA,lwd.tick=1),textOffsets=-.1)
  legend('bottom',sub('Other','Not VoC',sub('South Africa','SA',rownames(countTab))),fill=cols[1:nrow(countTab)],inset=-.33,ncol=6,xpd=NA,cex=.8)
dev.off()

