source('analyzeNewData.R')

options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)

mod <- rstan::stan_model("model.stan")
if(!exists('stan_sample'))stan_sample<-runStan(countTab[,,'random'],mod)
pdf('bayesTestNewData.pdf',width=6,height=6)
  plotStan(stan_sample,countTab[,,'random'],baseDate,cols)
dev.off()


modFlat <- rstan::stan_model("model_source_flat.stan")
if(!exists('stan_sample3'))stan_sample3<-runStan2(countTab[,,'random'],countTab[,,'s drop'],countTab[,,'vaccine'],modFlat,30000)

mod2 <- rstan::stan_model("model_source.stan")
#runStan2<-function(countTab,hospitalTab,dropTab,vaccineTab,mod,iter=2000){
if(!exists('stan_sample2'))stan_sample2<-runStan2(countTab[,,'random'],countTab[,,'s drop'],countTab[,,'vaccine'],mod2,30000)
meanLowUp<-calcCredInt(stan_sample2,names=rownames(countTab))
dense<-calcDensity(stan_sample2,names=rownames(countTab))

pdf('bayesTestData.pdf',width=7,height=6)
  layout(matrix(c(1,0,2,0,3,0),ncol=1),height=c(1,.1,1,.1,1,.2))
  par(mar=c(3,3,2,2),tcl=-.2,mgp=c(2,.3,0))
  plotBars(countTab[,,'random'],baseDate,cols,showLegend=FALSE)
  title(main='"Random" population',adj=.45,line=1)
  plotBars(countTab[,,'s drop'],baseDate,cols,showLegend=FALSE)
  title(main='S dropout',adj=.45,line=1)
  plotBars(countTab[,,'vaccine'],baseDate,cols,showLegend=FALSE)
  title(main='Vaccine breakthrough',adj=.45,line=1)
  legend('bottomleft',rownames(countTab),fill=cols[rownames(countTab)],inset=c(.1,-.6),ncol=4,xpd=NA,cex=.8)
dev.off()


pdf('bayesTestDataAndFit.pdf',width=7,height=6)
  layout(matrix(c(1,0,2,0),ncol=1),height=c(1,.1,1,.27))
  par(mar=c(3,3,2,2),tcl=-.2,mgp=c(2,.3,0))
  plotBars(meanLowUp$mean,baseDate,cols,showLegend=FALSE,showSum=FALSE)
  title(main='Predicted proportions data',adj=.45)
  plotBars(countTab[,,'random'],baseDate,cols,showLegend=FALSE)
  title(main='"Random" population',adj=.45,line=1)
  legend('bottomleft',rownames(countTab),fill=cols[rownames(countTab)],inset=c(.1,-.55),ncol=4,xpd=NA,cex=.8)
dev.off()

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

pdf('bayesTestProportionsMulti.pdf',width=7,height=4)
  par(mar=c(3,3,.5,3),tcl=-.2,mgp=c(2,.3,0))
  plotMeanUpperLower(meanLowUp[[1]],meanLowUp[[2]],meanLowUp[[3]],countTab[,,'random'],baseDate,cols)
dev.off()



mat<-as.matrix(stan_sample2)
zzz<-t(apply(mat[,grepl('dropChange',colnames(mat))],2,function(xx)c(mean(xx),quantile(xx,c(.025,.975)))))
print(round(t(t(apply(mat[,grepl('dropChange',colnames(mat))],2,function(xx)c(mean(xx<0))))),2),row.names=FALSE)
print(round(t(t(apply(mat[,grepl('vaccineChange',colnames(mat))],2,function(xx)c(mean(xx<0))))),2),row.names=FALSE)
rownames(zzz)<-rownames(countTab)
round(exp(zzz),1)


