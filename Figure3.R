
#Loading in packages
library(nlme)
library(mvtnorm)
library(ggplot2)
#library(lsmeans)
library(ggpubr)


######################################################
#Initializing functions.

#Creating correlation functions to generate residual correlation matrices for simulations.

cs.p<-function(rho,time.pts){
  mat.init<-matrix(rep(rho,times=time.pts*time.pts),time.pts,time.pts)
  diag.init<-diag(1-rho,time.pts)
  mat.final<-diag.init+mat.init
  return(mat.final)
}

ar.p<-function(rho,time.pts){
  t <- time.pts
  mat.final <- rho^abs(outer(t,t,"-"))
  return(mat.final)
}


gaus.p<-function(rho,time.pts){
  t <- time.pts
  mat.final <- rho^((outer(t,t,"-"))^2)
  return(mat.final)
}


#########################################################################
#Initializing simulation parameters

#The following code directly reads in the variability estimates from the GitHub website and
#stores them. User may also download the csv file directly from the GitHub page and read in the
#file from a local directory.
var.est<-read.csv("https://raw.githubusercontent.com/brakefieb/MiCor-LGE/main/varest.csv",stringsAsFactors = TRUE)
#var.est<-read.csv("/Users/jturner/Desktop/LongitudinalProject/Codes/varest.csv",stringsAsFactors = TRUE)
var.est<-var.est[,1]

#For curious reviewers
#summary(sqrt(var.est))
#hist(sqrt(var.est),breaks=100,xlim=c(0,2))


#Time point sampling for the design
n.time<-c(0,1,3,7,28)

#Number of subjects repeatedly sampled over time
n.sub<-10
n.sims=50 #don't change
rho<-0.2     
n.genes=10000 


#Creating study design of the repeated measures data set.
#Simulated gene expression values will be appended to this file for fitting linear models during
#the simulation loops.
des<-data.frame(time=rep(n.time,n.sub),sub=paste("S",rep(1:n.sub,each=length(n.time)),sep=""))
des$time2<-des$time
des$time<-factor(paste("T",des$time,sep=""))




###########################################################################
#Performing 3 separate simulations for each of the 3 scenarios in Figure 3

#Fig 3A  FDR Study (Truth is Compound Symmetry (CS))


#initializing result objects

disc<-rep(0,(length(n.time)-1)*3)  #Count vector for the total number of discoveries found after FDR correction
true.disc<-rep(0,(length(n.time)-1)*3)  #Count vector for the total number of true discoveries found after FDR correction


#On rare occasions, the gls model (using REML) will fail to converge, typically when fitting AR 
#or GAUS structures. We simply note them with a warning when they occur and move forward in the
#simulation loop.
set.seed(1234)
for(j in 1:n.sims){
  
  result.cs <- matrix(rep(0,10000*4),10000,4)
  result.sp <- matrix(rep(0,10000*4),10000,4)
  result.gaus<-matrix(rep(0,10000*4),10000,4)
  
  for (i in 1:n.genes){
    
    
    cor.mat<-cs.p(rho,time.pts=length(n.time))
    var.mat<-var.est[i]*cor.mat
    if(i < 201){resp.matrix<-rmvnorm(n.sub,mean=c(0,rep(1,length(n.time)-1)),sigma=var.mat)}
    if(i >200){resp.matrix<-rmvnorm(n.sub,mean=c(rep(0,length(n.time))),sigma=var.mat)}
    des$y<-as.vector(t(resp.matrix))
    
    
    tryCatch({
      test<-gls(model = y~ time, data = des, correlation = corCompSymm(form = ~1|sub),control = glsControl(opt = "optim")) 
      test2<-gls(model = y~ time, data = des, correlation = corExp(form = ~time2|sub),control = glsControl(opt = "optim"))
      test3<-gls(model = y~ time, data = des, correlation = corGaus(form = ~time2|sub),control = glsControl(opt = "optim"))
      
    },error=function(e){cat("Warning: Row",i,"\n")})
    
    
    result.cs[i,]<-coef(summary(test))[-1,4]
    result.sp[i,]<-coef(summary(test2))[-1,4]
    result.gaus[i,]<-coef(summary(test3))[-1,4]
    
  }
  
  adjp.cs<-apply(result.cs,2,p.adjust,method="fdr")
  adjp.sp<-apply(result.sp,2,p.adjust,method="fdr")
  adjp.gaus<-apply(result.gaus,2,p.adjust,method="fdr")
  
  disc.temp<-c(apply(adjp.cs,2,function(x){sum(x<.1)}),
               apply(adjp.sp,2,function(x){sum(x<.1)}),
               apply(adjp.gaus,2,function(x){sum(x<.1)}))
  
  true.disc.temp<-c(apply(adjp.cs[1:200,],2,function(x){sum(x<.1)}),
                    apply(adjp.sp[1:200,],2,function(x){sum(x<.1)}),
                    apply(adjp.gaus[1:200,],2,function(x){sum(x<.1)}))
  disc<-disc+disc.temp
  true.disc<-true.disc+true.disc.temp
}

#Calculation False Discovery and True Positive Rates
fdr<-(disc-true.disc)/disc
fdrME<-1.96*sqrt(fdr*(1-fdr)/disc) 

tpr<-true.disc/n.sims/200
tprME<-1.96*sqrt(tpr*(1-tpr)/(n.sims*200))

result<-data.frame(Struc=rep(c("CS","AR","Gaus"),each=4),FDR=fdr,ME=fdrME,TPR=tpr,ME_TPR=tprME)
result$time<- rep(c("T1 vs T0","T28 vs T0", "T3 vs T0", "T7 vs T0"),3)
result$time<-factor(result$time,levels=c("T1 vs T0", "T3 vs T0", "T7 vs T0","T28 vs T0"))


#write.csv(result,"FDR_CS_Rho2.csv")


fdr1<-ggplot(data=result, aes(x=time, y=FDR, fill=Struc))+geom_bar(stat="identity",position=position_dodge())+
  ylim(0,.45)+geom_errorbar(aes(ymin=FDR-ME,ymax=FDR+ME),width=.2,position=position_dodge(.9))+
  scale_fill_manual(values=c('coral1','deepskyblue1','chartreuse3'))+
  geom_hline(yintercept=0.1,size=0.75,linetype="longdash")+
  xlab("Comparison")+ylab("False Discovery Rate")+
  labs(fill='Structure') +ggtitle( bquote(CS~rho==0.2))+
  theme(legend.position="bottom",legend.background = element_rect(color = "black", 
                                                                  fill = "grey90", size = .75, linetype = "solid"),plot.title = element_text(hjust = 0.5))

fdr1

tpr1<-ggplot(data=result, aes(x=time, y=TPR, fill=Struc))+geom_bar(stat="identity",position=position_dodge())+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  geom_errorbar(aes(ymin=TPR-ME_TPR,ymax=TPR+ME_TPR),width=.2,position=position_dodge(.9))+
  scale_fill_manual(values=c('coral1','deepskyblue1','chartreuse3'))+
  geom_hline(yintercept=1,size=1.2,linetype="longdash")+
  xlab("Comparison")+ylab("True Positive Rate")+
  labs(fill='Structure') +ggtitle( bquote(CS~rho==0.2))+
  theme(legend.position="bottom",legend.background = element_rect(color = "black", 
                                                                  fill = "grey90", size = .75, linetype = "solid"),plot.title = element_text(hjust = 0.5))

tpr1








#Fig 3B  FDR Study (Truth is Autoregressivy (AR))


#initializing result objects

disc2<-rep(0,(length(n.time)-1)*3)  #Count vector for the total number of discoveries found after FDR correction
true.disc2<-rep(0,(length(n.time)-1)*3)  #Count vector for the total number of true discoveries found after FDR correction


#On rare occasions, the gls model (using REML) will fail to converge, typically when fitting AR 
#or GAUS structures. We simply note them with a warning when they occur and move forward in the
#simulation loop.
set.seed(1234)
for(j in 1:n.sims){
  
  result.cs <- matrix(rep(0,10000*4),10000,4)
  result.sp <- matrix(rep(0,10000*4),10000,4)
  result.gaus<-matrix(rep(0,10000*4),10000,4)
  
  for (i in 1:n.genes){
    
    cor.mat<-ar.p(rho,time.pts=n.time)
    var.mat<-var.est[i]*cor.mat
    if(i < 201){resp.matrix<-rmvnorm(n.sub,mean=c(0,rep(1,length(n.time)-1)),sigma=var.mat)}
    if(i >200){resp.matrix<-rmvnorm(n.sub,mean=c(rep(0,length(n.time))),sigma=var.mat)}
    des$y<-as.vector(t(resp.matrix))
    
    
    tryCatch({
      test<-gls(model = y~ time, data = des, correlation = corCompSymm(form = ~1|sub),control = glsControl(opt = "optim")) 
      test2<-gls(model = y~ time, data = des, correlation = corExp(form = ~time2|sub),control = glsControl(opt = "optim"))
      test3<-gls(model = y~ time, data = des, correlation = corGaus(form = ~time2|sub),control = glsControl(opt = "optim"))
      
    },error=function(e){cat("Warning: Row",i,"\n")})
    
    
    result.cs[i,]<-coef(summary(test))[-1,4]
    result.sp[i,]<-coef(summary(test2))[-1,4]
    result.gaus[i,]<-coef(summary(test3))[-1,4]
    
  }
  
  adjp.cs<-apply(result.cs,2,p.adjust,method="fdr")
  adjp.sp<-apply(result.sp,2,p.adjust,method="fdr")
  adjp.gaus<-apply(result.gaus,2,p.adjust,method="fdr")
  
  disc.temp<-c(apply(adjp.cs,2,function(x){sum(x<.1)}),
               apply(adjp.sp,2,function(x){sum(x<.1)}),
               apply(adjp.gaus,2,function(x){sum(x<.1)}))
  
  true.disc.temp<-c(apply(adjp.cs[1:200,],2,function(x){sum(x<.1)}),
                    apply(adjp.sp[1:200,],2,function(x){sum(x<.1)}),
                    apply(adjp.gaus[1:200,],2,function(x){sum(x<.1)}))
  disc2<-disc2+disc.temp
  true.disc2<-true.disc2+true.disc.temp
}

#Calculating FDR and TPR
fdr<-(disc2-true.disc2)/disc2
fdrME<-1.96*sqrt(fdr*(1-fdr)/disc2) 

tpr<-true.disc2/n.sims/200
tprME<-1.96*sqrt(tpr*(1-tpr)/(n.sims*200))

result2<-data.frame(Struc=rep(c("CS","AR","Gaus"),each=4),FDR=fdr,ME=fdrME,TPR=tpr,ME_TPR=tprME)
result2$time<- rep(c("T1 vs T0","T28 vs T0", "T3 vs T0", "T7 vs T0"),3)
result2$time<-factor(result$time,levels=c("T1 vs T0", "T3 vs T0", "T7 vs T0","T28 vs T0"))

#write.csv(result2,"FDR_AR_Rho2.csv")

fdr2<-ggplot(data=result2, aes(x=time, y=FDR, fill=Struc))+geom_bar(stat="identity",position=position_dodge())+
  ylim(0,.45)+geom_errorbar(aes(ymin=FDR-ME,ymax=FDR+ME),width=.2,position=position_dodge(.9))+
  scale_fill_manual(values=c('coral1','deepskyblue1','chartreuse3'))+
  geom_hline(yintercept=0.1,size=0.75,linetype="longdash")+
  xlab("Comparison")+ylab("False Discovery Rate")+
  labs(fill='Structure') +ggtitle( bquote(AR~rho==0.2))+
  theme(legend.position="bottom",legend.background = element_rect(color = "black", 
                                                                  fill = "grey90", size = .75, linetype = "solid"),plot.title = element_text(hjust = 0.5))

fdr2

tpr2<-ggplot(data=result2, aes(x=time, y=TPR, fill=Struc))+geom_bar(stat="identity",position=position_dodge())+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  geom_errorbar(aes(ymin=TPR-ME_TPR,ymax=TPR+ME_TPR),width=.2,position=position_dodge(.9))+
  scale_fill_manual(values=c('coral1','deepskyblue1','chartreuse3'))+
  geom_hline(yintercept=1,size=1.2,linetype="longdash")+
  xlab("Comparison")+ylab("True Positive Rate")+
  labs(fill='Structure') +ggtitle( bquote(AR~rho==0.2))+
  theme(legend.position="bottom",legend.background = element_rect(color = "black", 
                                                                  fill = "grey90", size = .75, linetype = "solid"),plot.title = element_text(hjust = 0.5))

tpr2





#Fig 3C  FDR Study (Truth is Autoregressivy (AR))


#initializing result objects

disc3<-rep(0,(length(n.time)-1)*3)  #Count vector for the total number of discoveries found after FDR correction
true.disc3<-rep(0,(length(n.time)-1)*3)  #Count vector for the total number of true discoveries found after FDR correction


#On rare occasions, the gls model (using REML) will fail to converge, typically when fitting AR 
#or GAUS structures. We simply note them with a warning when they occur and move forward in the
#simulation loop.
set.seed(1234)
for(j in 1:n.sims){
  
  result.cs <- matrix(rep(0,10000*4),10000,4)
  result.sp <- matrix(rep(0,10000*4),10000,4)
  result.gaus<-matrix(rep(0,10000*4),10000,4)
  
  for (i in 1:n.genes){
    
    cor.mat<-gaus.p(rho,time.pts=n.time)
    var.mat<-var.est[i]*cor.mat
    if(i < 201){resp.matrix<-rmvnorm(n.sub,mean=c(0,rep(1,length(n.time)-1)),sigma=var.mat)}
    if(i >200){resp.matrix<-rmvnorm(n.sub,mean=c(rep(0,length(n.time))),sigma=var.mat)}
    des$y<-as.vector(t(resp.matrix))
    
    
    tryCatch({
      test<-gls(model = y~ time, data = des, correlation = corCompSymm(form = ~1|sub),control = glsControl(opt = "optim")) 
      test2<-gls(model = y~ time, data = des, correlation = corExp(form = ~time2|sub),control = glsControl(opt = "optim"))
      test3<-gls(model = y~ time, data = des, correlation = corGaus(form = ~time2|sub),control = glsControl(opt = "optim"))
      
    },error=function(e){cat("Warning: Row",i,"\n")})
    
    
    result.cs[i,]<-coef(summary(test))[-1,4]
    result.sp[i,]<-coef(summary(test2))[-1,4]
    result.gaus[i,]<-coef(summary(test3))[-1,4]
    
  }
  
  adjp.cs<-apply(result.cs,2,p.adjust,method="fdr")
  adjp.sp<-apply(result.sp,2,p.adjust,method="fdr")
  adjp.gaus<-apply(result.gaus,2,p.adjust,method="fdr")
  
  disc.temp<-c(apply(adjp.cs,2,function(x){sum(x<.1)}),
               apply(adjp.sp,2,function(x){sum(x<.1)}),
               apply(adjp.gaus,2,function(x){sum(x<.1)}))
  
  true.disc.temp<-c(apply(adjp.cs[1:200,],2,function(x){sum(x<.1)}),
                    apply(adjp.sp[1:200,],2,function(x){sum(x<.1)}),
                    apply(adjp.gaus[1:200,],2,function(x){sum(x<.1)}))
  disc3<-disc3+disc.temp
  true.disc3<-true.disc3+true.disc.temp
}


fdr<-(disc3-true.disc3)/disc3
fdrME<-1.96*sqrt(fdr*(1-fdr)/disc3) 

tpr<-true.disc3/n.sims/200
tprME<-1.96*sqrt(tpr*(1-tpr)/(n.sims*200))

result3<-data.frame(Struc=rep(c("CS","AR","Gaus"),each=4),FDR=fdr,ME=fdrME,TPR=tpr,ME_TPR=tprME)
result3$time<- rep(c("T1 vs T0","T28 vs T0", "T3 vs T0", "T7 vs T0"),3)
result3$time<-factor(result$time,levels=c("T1 vs T0", "T3 vs T0", "T7 vs T0","T28 vs T0"))

#write.csv(result3,"FDR_GAUS_Rho2.csv")

fdr3<-ggplot(data=result3, aes(x=time, y=FDR, fill=Struc))+geom_bar(stat="identity",position=position_dodge())+
  ylim(0,.45)+geom_errorbar(aes(ymin=FDR-ME,ymax=FDR+ME),width=.2,position=position_dodge(.9))+
  scale_fill_manual(values=c('coral1','deepskyblue1','chartreuse3'))+
  geom_hline(yintercept=0.1,size=0.75,linetype="longdash")+
  xlab("Comparison")+ylab("False Discovery Rate")+
  labs(fill='Structure') +ggtitle( bquote(GAUS~rho==0.2))+
  theme(legend.position="bottom",legend.background = element_rect(color = "black", 
                                                                  fill = "grey90", size = .75, linetype = "solid"),plot.title = element_text(hjust = 0.5))

fdr3

tpr3<-ggplot(data=result3, aes(x=time, y=TPR, fill=Struc))+geom_bar(stat="identity",position=position_dodge())+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  geom_errorbar(aes(ymin=TPR-ME_TPR,ymax=TPR+ME_TPR),width=.2,position=position_dodge(.9))+
  scale_fill_manual(values=c('coral1','deepskyblue1','chartreuse3'))+
  geom_hline(yintercept=1,size=1.2,linetype="longdash")+
  xlab("Comparison")+ylab("True Positive Rate")+
  labs(fill='Structure') +ggtitle( bquote(GAUS~rho==0.2))+
  theme(legend.position="bottom",legend.background = element_rect(color = "black", 
                                                                  fill = "grey90", size = .75, linetype = "solid"),plot.title = element_text(hjust = 0.5))

tpr3




#Producing final figure
#png("Figure3.png", width = 10, height = 4.25, units="in", res=1200)
ggarrange(fdr1,fdr2,fdr3,ncol=3,nrow=1,common.legend=T, labels=c("A","B","C"),legend="bottom")
#dev.off()


#Producing final figure
#png("Supp2.png", width = 10, height = 4.25, units="in", res=1200)
ggarrange(tpr1,tpr2,tpr3,ncol=3,nrow=1,common.legend=T, labels=c("A","B","C"),legend="bottom")
#dev.off()
