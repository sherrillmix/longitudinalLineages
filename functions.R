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

calcCredInt<-function(stan_sample,quants=c(.025,.975),names=NULL,colPattern='^props\\['){
  mat<-as.matrix(stan_sample)
  props<-mat[,grepl(colPattern,colnames(mat))]
  meanLowerUpper<-apply(props,2,meanCrI,quants)
  convertToMat<-function(xx){
    row<-as.numeric(sub('[^0-9]+([0-9]+),[0-9]+.*','\\1',names(xx)))
    col<-as.numeric(sub('[^0-9]+[0-9]+,([0-9]+).*','\\1',names(xx)))
    tapply(xx,list(row,col),c)
  }
  means<-convertToMat(meanLowerUpper[1,])
  lower<-convertToMat(meanLowerUpper[2,])
  upper<-convertToMat(meanLowerUpper[3,])
  rownames(means)<-rownames(lower)<-rownames(upper)<-names
  return(list('mean'=means,'upper'=upper,'lower'=lower))
}

calcDensity<-function(stan_sample,names=NULL){
  mat<-as.matrix(stan_sample)
  drops<-mat[,grepl('^dropChange\\[',colnames(mat))]
  vaccines<-mat[,grepl('^vaccineChange\\[',colnames(mat))]
  minMax<-range(cbind(drops,vaccines))
  vaccineDensity<-apply(vaccines,2,density,from=minMax[1],to=minMax[2])
  dropDensity<-apply(drops,2,density,from=minMax[1],to=minMax[2])
  names(vaccineDensity)<-names(dropDensity)<-names
  return(list('vaccine'=vaccineDensity,'drop'=dropDensity))
}

plotMeanUpperLower<-function(means,upper,lower,countTab,baseDate,cols=NULL,...){
  #https://sashamaps.net/docs/resources/20-colors/
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  cols<-cols[names(cols) %in% rownames(countTab)]
  plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Estimated proportion',xaxs='i',...)
  for(ii in 1:nrow(means)){
    polygon(c(1:ncol(means),ncol(means):1),c(lower[ii,],rev(upper[ii,])),border=NA,col=sprintf('%s33',cols[rownames(countTab)[ii]]))
  }
  for(ii in 1:nrow(means)){
    lines(1:ncol(means),means[ii,],col=cols[rownames(countTab)[ii]])
  }
  abline(h=0)
  text(ncol(means)+.5,means[,ncol(means)],rownames(countTab),cex=.7,xpd=NA,adj=0)
  prettyDates<-lubridate::mdy(c('4/1/20','7/1/20','10/1/20','1/1/21','4/1/21'))
  axis(1,(prettyDates-baseDate)/7+1,sub('  +',' ',format(prettyDates,'%b %e %Y')),cex=.8,tcl=-.5,mgp=c(3,.7,0))
}

plotBars<-function(countTab,baseDate,cols=NULL,showLegend=TRUE,showSum=TRUE,showDates=TRUE,skipDate=NULL,dateCex=.6){
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  cols<-cols[names(cols) %in% rownames(countTab)]
  counts<-apply(countTab,2,sum)
  propTab<-apply(countTab,2,function(xx)xx/sum(xx))
  plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-1,1),ylim=c(0,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Proportion',xaxs='i')
  for(ii in 1:ncol(propTab))rect(ii-.5,cumsum(propTab[,ii]),ii+.5,cumsum(c(0,propTab[-nrow(propTab),ii])),col=cols[rownames(propTab)],lwd=.5)
  uniqDate<-data.frame('week'=which(counts>0),'rdate'=baseDate+(which(counts>0)-1)*7)
  if(is.null(skipDate))dnar::slantAxis(1,uniqDate$week,if(showDates)sub('  +',' ',format(uniqDate$rdate,'%b %e %Y')) else rep('',nrow(uniqDate)),location=.5,cex=dateCex)
  else dnar::slantAxis(1,uniqDate$week,ifelse(1:nrow(uniqDate)%%skipDate==1,sub('  +',' ',format(uniqDate$rdate,'%b %e %Y')),''),location=.5,cex=dateCex)
  if(showSum)dnar::slantAxis(3,which(counts>0),counts[counts>0],location=.4,cex=.6,axisArgs=list(lwd=NA,lwd.tick=1),textOffsets=-.1)
  if(showLegend)legend('bottom',names(cols),fill=cols,inset=-.33,ncol=ceiling(nrow(countTab)/2),xpd=NA,cex=.8)
}


runStan2<-function(countTab,dropTab,vaccineTab,mod,iter=2000){#,hospitalTab
  #should do more double checking
  stan_sample <- rstan::sampling(
    mod,
    data=list(
      counts=countTab,
      #hospital=hospitalTab,
      drop=dropTab,
      vaccine=vaccineTab,
      nTime=ncol(countTab),
      nLineage=nrow(countTab)
    ),
    iter=iter,
    chains=50,
    thin=2,
    control=list(max_treedepth=15),
    pars=c('means','propsDrop','propsVaccine'),
    include=FALSE
  )
  return(stan_sample)
}

plotIndivStan<-function(means,upper,lower,countTab,baseDate,cols=NULL,nCol=5,additional=NULL){
  propTab<-apply(countTab,2,function(xx)xx/sum(xx))
  counts<-apply(countTab,2,sum)
  barCol<-rev(grey.colors(max(counts)+1))
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  #https://sashamaps.net/docs/resources/20-colors/
  prettyDates<-lubridate::mdy(c('4/1/20','7/1/20','10/1/20','1/1/21','4/1/21','7/1/21'))
  for(ii in 1:nrow(means)){
    plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-4,4),ylim=c(0,1.2),las=1,xlab='',bty='l',xaxt='n',yaxs='i',xaxs='i',yaxt='n',xaxt='n')
    prettyY<-c(0,.5,1)
    if(ii %% nCol==1)axis(2,las=1,prettyY)
    else axis(2,prettyY,rep('',length(prettyY)))
    if(ii > nrow(countTab)-nCol) dnar::slantAxis(1,(prettyDates-baseDate)/7+1,sub('  +',' ',format(prettyDates,'%b %e %Y')),cex=.8,tcl=-.2,location=.5)
    else axis(1,(prettyDates-baseDate)/7+1,rep('',length(prettyDates)),tcl=-.2)
    if(ii==(floor(nrow(means)/nCol/2)*nCol+1))mtext('Estimated proportion',2,line=2.2,at=1)
    title(main=rownames(countTab)[ii],line=-1.3,cex=.9)
    rect(1:ncol(propTab)-.5,0,1:ncol(propTab)+.5,propTab[ii,],col=barCol[counts+1],border=NA)
    polygon(c(1:ncol(means),ncol(means):1),c(lower[ii,],rev(upper[ii,])),border=NA,col=sprintf('%s77',cols[rownames(countTab)[ii]]))
    lines(1:ncol(means),means[ii,],col=cols[rownames(countTab)[ii]])
    if(!is.null(additional)){
      lines(1:ncol(additional$mean),additional$mean[ii,],col=cols[rownames(countTab)[ii]],lty=2)
      polygon(c(1:ncol(additional$mean),ncol(additional$mean):1),c(additional$lower[ii,],rev(additional$upper[ii,])),border=NA,col=sprintf('%s66',cols[rownames(countTab)[ii]]))
    }
    box()
  }
  abline(h=0)
}


plotIndivDense<-function(dense,cols=NULL,xlim=range(exp(dense[[1]]$x)),ylab='Fold change'){
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  #https://sashamaps.net/docs/resources/20-colors/
  yMax<-max(sapply(dense,function(xx)max(xx$y)))
  for(ii in 1:length(dense)){
    plot(1,1,type='n',xlim=xlim,ylim=c(0,yMax),las=1,xlab='',bty='n',xaxt='n',yaxs='i',xaxs='i',yaxt='n',xaxt='n',log='x')
    polygon(exp(dense[[ii]]$x),dense[[ii]]$y,col=cols[names(dense)[ii]])
    title(main=names(dense)[ii],line=0,cex=.9)
    dnar::logAxis(1,axisVals=c(-2,0,2,4))
    abline(v=1,lty=2)
    #if(ii==(2*nCol+1))mtext(2,line=2.2,at=par('usr')[4])
  }
  text(grconvertX(.015,'ndc','user'),grconvertY(.5,'ndc','user'),'Estimated posterior probability',xpd=NA,cex=1.2,srt=90)
  text(grconvertX(.5,'ndc','user'),grconvertY(.02,'ndc','user'),ylab,xpd=NA,cex=1.2)
}

#deprecated
plotStan<-function(stan_sample,countTab,baseDate,cols=NULL){
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  cols<-cols[names(cols) %in% rownames(countTab)]
  counts<-apply(countTab,2,sum)
  propTab<-apply(countTab,2,function(xx)xx/sum(xx))
  mat<-as.matrix(stan_sample)
  props<-mat[,grepl('^props\\[',colnames(mat))]
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
    polygon(c(1:ncol(means),ncol(means):1),c(lower[ii,],rev(upper[ii,])),border=NA,col=sprintf('%s33',cols[rownames(countTab)[ii]]))
  }
  for(ii in 1:nrow(means)){
    lines(1:ncol(means),means[ii,],col=cols[rownames(countTab)[ii]])
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

runStan3<-function(countTab,dropTab,vaccineTab,greek,mod,iter=2000){
  #should do more double checking
  if(length(greek)!=nrow(countTab))stop('Sublineage groupings not same length as count table')
  abundantGreek<-names(table(greek)[table(greek)>1])
  greekIds<-structure(1:(length(abundantGreek)+1),.Names=c('__BASE__',abundantGreek))
  dat<-list(
    counts=countTab,
    drop=dropTab,
    vaccine=vaccineTab,
    nTime=ncol(countTab),
    nLineage=nrow(countTab),
    nSublineageGroup=max(greekIds),
    sublineageGroup=greekIds[ifelse(greek %in% names(greekIds),greek,'__BASE__')]
  )
  stan_sample <- rstan::sampling(
    mod,
    data=dat,
    iter=iter,
    chains=50,
    thin=2,
    control=list(max_treedepth=15),
    pars=c('means'),
    include=FALSE
  )
  return(list('stan'=stan_sample,'dat'=dat,'greek'=greekIds,'lineage'=rownames(countTab)))
}

calcDensity2<-function(stan_sample,names=NULL,greekNames=NULL){
  mat<-as.matrix(stan_sample)
  drops<-mat[,grepl('^dropChange\\[',colnames(mat))]
  vaccines<-mat[,grepl('^vaccineChangeWithSub\\[',colnames(mat))]
  greeks<-mat[,grepl('^vaccineChange\\[',colnames(mat))]
  minMax<-range(cbind(drops,vaccines,greeks))
  vaccineDensity<-apply(vaccines,2,density,from=minMax[1],to=minMax[2])
  dropDensity<-apply(drops,2,density,from=minMax[1],to=minMax[2])
  greekDensity<-apply(greeks,2,density,from=minMax[1],to=minMax[2])
  names(vaccineDensity)<-names(dropDensity)<-names
  names(greekDensity)<-greekNames
  return(list('vaccine'=vaccineDensity,'drop'=dropDensity,'greek'=greekDensity))
}

runMutationStan<-function(countTab,dropTab,vaccineTab,mod,nChain=50,nIter=2000){
  dat<-list(
    counts=countTab[1,],
    nCounts=apply(countTab,2,sum),
    drop=dropTab[1,],
    nDrop=apply(dropTab,2,sum),
    vaccine=vaccineTab[1,],
    nVaccine=apply(vaccineTab,2,sum),
    nTime=ncol(countTab)
  )
  stan_sample <- rstan::sampling(
    mod,
    data=dat,
    iter=nIter,
    chains=nChain,
    thin=2,
    control=list(max_treedepth=15) #adapt_delta=.99
    #pars=c('means'),
    #include=FALSE
  )
  return(list('stan'=stan_sample,'dat'=dat,'mod'=mod))
}

logit<-function(xx)log(xx)-log(1-xx)
invLogit<-function(xx)1/(1+exp(-xx))
meanCrI<-function(xx,quants=c(.025,.975))c(mean(xx),quantile(xx,quants))


getTable<-function(stan,addGroup=TRUE){
  mat<-as.matrix(stan$stan)
  strainChange<-exp(apply(mat[,grepl('vaccineChangeWithSub\\[',colnames(mat))],2,meanCrI))
  strainChange<-rbind(strainChange,'prob'=apply(mat[,grepl('vaccineChangeWithSub\\[',colnames(mat))]>0,2,mean))
  colnames(strainChange)<-stan$lineage
  groupChange<-exp(apply(mat[,grepl('vaccineChangeGroup\\[',colnames(mat))],2,meanCrI))
  groupChange<-rbind(groupChange,'prob'=apply(mat[,grepl('vaccineChangeGroup\\[',colnames(mat))]>0,2,mean))
  colnames(groupChange)<-sub('Other','Non-VoI/VoC',names(stan$greek)[-1])
  if(addGroup)out<-t(cbind(groupChange,strainChange))
  else out<-t(strainChange)
  return(out)
}

getTableDrop<-function(stan){
  mat<-as.matrix(stan$stan)
  strainChange<-exp(apply(mat[,grepl('dropChange\\[',colnames(mat))],2,meanCrI))
  strainChange<-rbind(strainChange,'prob'=apply(mat[,grepl('dropChange\\[',colnames(mat))]>0,2,mean))
  colnames(strainChange)<-stan$lineage
  out<-t(strainChange)
  return(out)
}

plotMeanUpperLowerLog<-function(means,upper,lower,countTab,baseDate,cols=NULL,...){
  #https://sashamaps.net/docs/resources/20-colors/
  if(is.null(cols)) cols<-structure(c('#469990','#e6194B', '#4363d8', '#ffe119', '#f58231', '#42d4f4', '#3cb44b','#469990','#800000','#808000','#911eb4')[1:nrow(countTab)],.Names=rownames(countTab))
  cols<-cols[names(cols) %in% rownames(countTab)]
  for(ii in 1:nrow(means)){
    plot(1,1,type='n',xlim=c(1,ncol(countTab))+c(-1,1),ylim=c(.001,1),las=1,xlab='',bty='l',xaxt='n',yaxs='i',ylab='Estimated proportion',xaxs='i',log='y',main=rownames(countTab)[ii],...)
    polygon(c(1:ncol(means),ncol(means):1),c(lower[ii,],rev(upper[ii,])),border=NA,col=sprintf('%s33',cols[rownames(countTab)[ii]]))
    lines(1:ncol(means),means[ii,],col=cols[rownames(countTab)[ii]])
    prettyDates<-lubridate::mdy(c('4/1/20','7/1/20','10/1/20','1/1/21','4/1/21'))
    axis(1,(prettyDates-baseDate)/7+1,sub('  +',' ',format(prettyDates,'%b %e %Y')),cex=.8,tcl=-.5,mgp=c(3,.7,0))
  }
}
