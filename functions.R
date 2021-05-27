runStan<-function(countTab,mod,iter=2000){
  stan_sample <- rstan::sampling(
    mod,
    data=list(
      counts=countTab,
      nTime=ncol(countTab),
      nLineage=nrow(countTab)
    ),
    iter=iter,
    chains=50,
    control=list(max_treedepth=15),
  )
  return(stan_sample)
}

plotStan<-function(stan_sample,countTab,baseDate,cols=NULL){
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  cols<-cols[names(cols) %in% rownames(countTab)]
  counts<-apply(countTab,2,sum)
  propTab<-apply(countTab,2,function(xx)xx/sum(xx))
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
  par(mar=c(3,3,.5,4),tcl=-.2,mgp=c(2,.3,0))
  layout(matrix(c(1,0,2,0),ncol=1),height=c(1,.01,1,.15))
  plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Estimated proportion',xaxs='i')
  for(ii in 1:nrow(means)){
    lines(1:ncol(means),means[ii,],col=cols[rownames(countTab)[ii]])
    polygon(c(1:ncol(means),ncol(means):1),c(lower[ii,],rev(upper[ii,])),border=NA,col=sprintf('%s33',cols[rownames(countTab)[ii]]))
  }
  abline(h=0)
  text(ncol(means)+.5,means[,ncol(means)],rownames(countTab),cex=.7,xpd=NA,adj=0)
  prettyDates<-lubridate::mdy(c('4/1/20','7/1/20','10/1/20','1/1/21','4/1/21'))
  axis(1,(prettyDates-baseDate)/7+1,sub('  +',' ',format(prettyDates,'%b %e %Y')),cex=.8,tcl=-.5,mgp=c(3,.7,0))
  #
  plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Proportion',xaxs='i')
  for(ii in 1:ncol(propTab))rect(ii-.5,cumsum(propTab[,ii]),ii+.5,cumsum(c(0,propTab[-nrow(propTab),ii])),col=cols[rownames(propTab)])
  uniqDate<-data.frame('week'=which(counts>0),'rdate'=baseDate+(which(counts>0)-1)*7)
  dnar::slantAxis(1,uniqDate$week,sub('  +',' ',format(uniqDate$rdate,'%b %e %Y')),location=.5,cex=.6)
  dnar::slantAxis(3,which(counts>0),counts[counts>0],location=.4,cex=.6,axisArgs=list(lwd=NA,lwd.tick=1),textOffsets=-.1)
  legend('bottom',names(cols),fill=cols,inset=-.33,ncol=ceiling(nrow(countTab)/2),xpd=NA,cex=.8)
  invisible(return(list(means,lower,upper)))
}

