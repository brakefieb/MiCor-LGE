


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
var.est<-read.csv("https://raw.githubusercontent.com/brakefieb/Misspecification-of-Correlation-Structures-in-Longitudinal-Gene-Expression-Studies-R-Scripts/main/varest.csv",stringsAsFactors = TRUE)
#var.est<-read.csv("/Users/jturner/Desktop/LongitudinalProject/Codes/varest.csv",stringsAsFactors = TRUE)
var.est<-var.est[,1]

#For curious reviewers
#summary(sqrt(var.est))
#hist(sqrt(var.est),breaks=100,xlim=c(0,2))


#Time point sampling for the design
n.time<-c(0,1,3,7,28)

#Number of subjects repeatedly sampled over time
n.sub<-10
n.sims=10000 #don't change
rho<-0.2     


#Creating study design of the repeated measures data set.
#Simulated gene expression values will be appended to this file for fitting linear models during
#the simulation loops.
des<-data.frame(time=rep(n.time,n.sub),sub=paste("S",rep(1:n.sub,each=length(n.time)),sep=""))
des$time2<-des$time
des$time<-factor(paste("T",des$time,sep=""))
View(des)






###########################################################################
#Performing 3 separate simulations for each of the 3 scenarios in Figure 2

#Fig 2A  Type-I Error Study (Truth is Compound Symmetry (CS))

#Initializing result vectors that record p-values for each of the comparisons and across all 10000 tests
result.cs <- c()
result.ar <- c()
result.gaus<-c()

#On rare occasions, the gls model (using REML) will fail to converge, typically when fitting AR 
#or GAUS structures. We simply note them with a warning when they occur and move forward in the
#simulation loop.
set.seed(1234)
for (i in 1:n.sims){
  
  cor.mat<-cs.p(rho,time.pts=length(n.time))                            #creating residual side correlation matrix for each subject
  var.mat<-var.est[i]*cor.mat                                           #creating residual side covariance matrix for each subject
  resp.matrix<-rmvnorm(n.sub,mean=rep(0,length(n.time)),sigma=var.mat)  #Generating gene epression values with specified correlation 
  des$y<-as.vector(t(resp.matrix))                                      #Appending response values to design file for modeling
  
  
  tryCatch({
    model.cs<-gls(model = y~ time, data = des, correlation = corCompSymm(form = ~1|sub),control = glsControl(opt = "optim")) 
    model.ar<-gls(model = y~ time, data = des, correlation = corExp(form = ~time2|sub),control = glsControl(opt = "optim"))
    model.gaus<-gls(model = y~ time, data = des, correlation = corGaus(form = ~time2|sub),control = glsControl(opt = "optim"))
    
  },error=function(e){cat("Warning: Row",i,"\n")})
  
  
  result.cs<-rbind(result.cs,coef(summary(model.cs))[-1,4])
  result.ar<-rbind(result.ar,coef(summary(model.ar))[-1,4])
  result.gaus<-rbind(result.gaus,coef(summary(model.gaus))[-1,4])
  
}


type1.cs<-apply(result.cs,2,function(x){sum(x<.05)/length(x)})
type1.ar<-apply(result.ar,2,function(x){sum(x<.05)/length(x)})
type1.gaus<-apply(result.gaus,2,function(x){sum(x<.05)/length(x)})



result1<-data.frame(Struc=rep(c("CS","AR","Gaus"),each=4),error=c(type1.cs,type1.ar,type1.gaus))
result1$time<- rep(c("T1 vs T0","T28 vs T0", "T3 vs T0", "T7 vs T0"),3)
result1$time<-factor(result1$time,levels=c("T1 vs T0", "T3 vs T0", "T7 vs T0","T28 vs T0"))
result1$ME<-1.96*sqrt(result1$error*(1-result1$error)/10000)



ggplot(data=result1, aes(x=time, y=error, fill=Struc))+geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=error-ME,ymax=error+ME),width=.2,position=position_dodge(.9))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values=c('coral1','deepskyblue1','chartreuse3'))+
  geom_hline(yintercept=0.05,size=0.75,linetype="longdash")+
  xlab("Comparison")+ylab("Type-I Error")+
  labs(fill='Structure') +ggtitle( bquote(CS~rho==0.2))+
  theme(legend.position="bottom",legend.background = element_rect(color = "black", 
                                                                  fill = "grey90", size = .75, linetype = "solid"),plot.title = element_text(hjust = 0.5))

#Storing plot to make final figure
p1<-ggplot(data=result1, aes(x=time, y=error, fill=Struc))+geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=error-ME,ymax=error+ME),width=.2,position=position_dodge(.9))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values=c('coral1','deepskyblue1','chartreuse3'))+
  geom_hline(yintercept=0.05,size=0.75,linetype="longdash")+
  xlab("Comparison")+ylab("Type-I Error")+
  labs(fill='Structure') +ggtitle( bquote(CS~rho==0.2))+
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))




#Fig 2B  Type-I Error Study (Truth is Autoregressive (AR))

#Initializing result vectors that record p-values for each of the comparisons and across all 10000 tests
result.cs2 <- c()
result.ar2 <- c()
result.gaus2 <- c()


#On rare occasions, the gls model (using REML) will fail to converge, typically when fitting AR 
#or GAUS structures. We simply note them with a warning when they occur and move forward in the
#simulation loop.
set.seed(1234)
for (i in 1:n.sims){
  cor.mat<-ar.p(rho,time.pts=n.time)
  var.mat<-var.est[i]*cor.mat
  resp.matrix<-rmvnorm(n.sub,mean=rep(0,length(n.time)),sigma=var.mat)
  des$y<-as.vector(t(resp.matrix))
  
  tryCatch({
    model.cs2<-gls(model = y~ time, data = des, correlation = corCompSymm(form = ~1|sub),control = glsControl(opt = "optim")) 
    model.ar2<-gls(model = y~ time, data = des, correlation = corExp(form = ~time2|sub),control = glsControl(opt = "optim"))
    model.gaus2<-gls(model = y~ time, data = des, correlation = corGaus(form = ~time2|sub),control = glsControl(opt = "optim"))
  },error=function(e){cat("Warning: Row",i,"\n")})
  
  
  result.cs2<-rbind(result.cs2,coef(summary(model.cs2))[-1,4])
  result.ar2<-rbind(result.ar2,coef(summary(model.ar2))[-1,4])
  result.gaus2<-rbind(result.gaus2,coef(summary(model.gaus2))[-1,4])
  
}

type1.cs2<-apply(result.cs2,2,function(x){sum(x<.05)/length(x)})
type1.ar2<-apply(result.ar2,2,function(x){sum(x<.05)/length(x)})
type1.gaus2<-apply(result.gaus2,2,function(x){sum(x<.05)/length(x)})

result2<-data.frame(Struc=rep(c("CS","AR","Gaus"),each=4),error=c(type1.cs2,type1.ar2,type1.gaus2))
result2$time<- rep(c("T1 vs T0","T28 vs T0", "T3 vs T0", "T7 vs T0"),3)
result2$time<-factor(result2$time,levels=c("T1 vs T0", "T3 vs T0", "T7 vs T0","T28 vs T0"))
result2$ME<-1.96*sqrt(result2$error*(1-result2$error)/10000)


ggplot(data=result2, aes(x=time, y=error, fill=Struc))+geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=error-ME,ymax=error+ME),width=.2,position=position_dodge(.9))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values=c('coral1','deepskyblue1','chartreuse3'))+
  geom_hline(yintercept=0.05,size=0.75,linetype="longdash")+
  xlab("Comparison")+ylab("Type-I Error")+
  labs(fill='Structure') +ggtitle( bquote(AR~rho==0.2))+
  theme(legend.position="bottom",legend.background = element_rect(color = "black", 
                                                                  fill = "grey90", size = .75, linetype = "solid"),plot.title = element_text(hjust = 0.5))

#Storing plot to make final figure
p2<-ggplot(data=result2, aes(x=time, y=error, fill=Struc))+geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=error-ME,ymax=error+ME),width=.2,position=position_dodge(.9))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values=c('coral1','deepskyblue1','chartreuse3'))+
  geom_hline(yintercept=0.05,size=0.75,linetype="longdash")+
  xlab("Comparison")+ylab("Type-I Error")+
  labs(fill='Structure') +ggtitle( bquote(AR~rho==0.2))+
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))



#Fig 2C  Type-I Error Study (Truth is Gaussian (GAUS))

#Initializing result vectors that record p-values for each of the comparisons and across all 10000 tests
result.cs3 <- c()
result.ar3 <- c()
result.gaus3 <- c()



#On rare occasions, the gls model (using REML) will fail to converge, typically when fitting AR 
#or GAUS structures. We simply note them with a warning when they occur and move forward in the
#simulation loop.
set.seed(1234)
for (i in 1:n.sims){
  cor.mat<-gaus.p(rho,time.pts=n.time)
  var.mat<-var.est[i]*cor.mat
  resp.matrix<-rmvnorm(n.sub,mean=rep(0,length(n.time)),sigma=var.mat)
  des$y<-as.vector(t(resp.matrix))
  
  tryCatch({
    model.cs3<-gls(model = y~ time, data = des, correlation = corCompSymm(form = ~1|sub),control = glsControl(opt = "optim")) 
    model.ar3<-gls(model = y~ time, data = des, correlation = corExp(form = ~time2|sub),control = glsControl(opt = "optim"))
    model.gaus3<-gls(model = y~ time, data = des, correlation = corGaus(form = ~time2|sub),control = glsControl(opt = "optim"))
  },error=function(e){cat("Warning: Row",i,"\n")})
  
  result.cs3<-rbind(result.cs3,coef(summary(model.cs3))[-1,4])
  result.ar3<-rbind(result.ar3,coef(summary(model.ar3))[-1,4])
  result.gaus3<-rbind(result.gaus3,coef(summary(model.gaus3))[-1,4])
  
}

type1.cs3<-apply(result.cs3,2,function(x){sum(x<.05)/length(x)})
type1.ar3<-apply(result.ar3,2,function(x){sum(x<.05)/length(x)})
type1.gaus3<-apply(result.gaus3,2,function(x){sum(x<.05)/length(x)})

result3<-data.frame(Struc=rep(c("CS","AR","Gaus"),each=4),error=c(type1.cs3,type1.ar3,type1.gaus3))
result3$time<- rep(c("T1 vs T0","T28 vs T0", "T3 vs T0", "T7 vs T0"),3)
result3$time<-factor(result3$time,levels=c("T1 vs T0", "T3 vs T0", "T7 vs T0","T28 vs T0"))
result3$ME<-1.96*sqrt(result3$error*(1-result3$error)/10000)


ggplot(data=result3, aes(x=time, y=error, fill=Struc))+geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=error-ME,ymax=error+ME),width=.2,position=position_dodge(.9))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values=c('coral1','deepskyblue1','chartreuse3'))+
  geom_hline(yintercept=0.05,size=0.75,linetype="longdash")+
  xlab("Comparison")+ylab("Type-I Error")+
  labs(fill='Structure') +ggtitle( bquote(Gaus~rho==0.2))+
  theme(legend.position="bottom",legend.background = element_rect(color = "black", 
                                                                  fill = "grey90", size = .75, linetype = "solid"),plot.title = element_text(hjust = 0.5))


p3<-ggplot(data=result3, aes(x=time, y=error, fill=Struc),palette="jco")+geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=error-ME,ymax=error+ME),width=.2,position=position_dodge(.9))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  scale_fill_manual(values=c('coral1','deepskyblue1','chartreuse3'))+
  geom_hline(yintercept=0.05,size=0.75,linetype="longdash")+
  xlab("Comparison")+ylab("Type-I Error")+
  labs(fill='Structure') +ggtitle( bquote(Gaus~rho==0.2))+
  theme(legend.position="none",plot.title = element_text(hjust = 0.5))


#Producing final figure
#png("Figure2.png", width = 10, height = 4.25, units="in", res=1200)
ggarrange(p1,p2,p3,ncol=3,nrow=1,common.legend=T, labels=c("A","B","C"),legend="bottom")
#dev.off()





































