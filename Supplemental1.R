library(nlme)
library(mvtnorm)
library(ggplot2)

######################################################
#Initializing functions.

#Creating correlation function to generate residual correlation matrices for an AR simulation.

ar.p<-function(rho,time.pts){
  t <- time.pts
  mat.final <- rho^abs(outer(t,t,"-"))
  return(mat.final)
}


#########################################################################
#Initializing simulation parameters

#The following code directly reads in the variability estimates from the GitHub website and
#stores them. User may also download the csv file directly from the GitHub page and read in the
#file from a local directory.
var.est<-read.csv("https://raw.githubusercontent.com/brakefieb/Misspecification-of-Correlation-Structures-in-Longitudinal-Gene-Expression-Studies-R-Scripts/main/varest.csv")
#var.est<-read.csv("/Users/jturner/Desktop/LongitudinalProject/Codes/varest.csv")
var.est<-var.est[,1]

#For curious reviewers
#summary(sqrt(var.est))
#hist(sqrt(var.est),breaks=100,xlim=c(0,2))


#Time point sampling for the design 
#(Adding additional equally spaced timepoints to illustrate that type-I error can be controlled
#when properly specified under the AR structure.)
n.time<-c(0,1,2,3,4,7,14,28)

#Number of subjects repeatedly sampled over time
n.sub<-10
n.sims=10000 #don't change
rho<-0.2     


#Creating study design of the repeated measures data set.
#Simulated gene expression values will be appended to this file for fitting linear models during
#the simulation loops.
des<-data.frame(time=rep(n.time,n.sub),sub=paste("S",rep(1:n.sub,each=length(n.time)),sep=""))
des$time2<-des$time
des$time<-factor(paste("T",des$time,sep=""),levels=c("T0","T1","T2","T3","T4","T7", "T14", "T28"))      


#Supplementary Figure 1:    Type-I Error Study (Truth is Autoregressive (AR) with more equally spaced time points)

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

result2<-data.frame(Struc=rep(c("CS","AR","Gaus"),each=7),error=c(type1.cs2,type1.ar2,type1.gaus2))
result2$time<- rep(c("T1 vs T0","T2 vs T0", "T3 vs T0","T4 vs T0", "T7 vs T0","T14 vs T0","T28 vs T0"),3)
result2$time<-factor(result2$time,levels=c("T1 vs T0","T2 vs T0", "T3 vs T0","T4 vs T0", "T7 vs T0","T14 vs T0","T28 vs T0"))
result2$ME<-1.96*sqrt(result2$error*(1-result2$error)/10000)
#library(latex2exp)

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
  theme(legend.position="bottom",legend.background = element_rect(color = "black", 
                                                                  fill = "grey90", size = .75, linetype = "solid"),plot.title = element_text(hjust = 0.5))



#png("Supp1.png", width = 6, height = 3.5, units="in", res=1200)
p2
#dev.off()


