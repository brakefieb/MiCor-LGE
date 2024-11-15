library(nlme)
library(mvtnorm)
library(ggplot2)
#library(lsmeans)
library(reshape2)
library(Hmisc)



#The following code directly reads in the variability estimates from the GitHub website and
#stores them. User may also download the csv file directly from the GitHub page and read in the
#file from a local directory.
var.est<-read.csv("https://raw.githubusercontent.com/brakefieb/MiCor-LGE/main/varest.csv",stringsAsFactors = TRUE)
#var.est<-read.csv("/Users/jturner/Desktop/LongitudinalProject/Codes/varest.csv",stringsAsFactors = TRUE)
var.est<-var.est[,1]


#Creating functions to generate covariance/correlation matrices for simulations.
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


#Initializing result objects
disc<-rep(0,(length(n.time)-1)*2)  #Count vector for the total number of discoveries found after FDR correction
true.disc<-rep(0,(length(n.time)-1)*2)  #Count vector for the total number of true discoveries found after FDR correction


#On rare occasions, the gls model (using REML) will fail to converge, typically when fitting AR 
#or GAUS structures. We simply note them with a warning when they occur and move forward in the
#simulation loop.
set.seed(1234)

for(k in 1:n.sims){
  #initializing result objects
  result.ar <- c()
  result.arp <- c()
  #Will store Expression Values in a matrix to compute emperical correlation function estimate
  eset<-matrix(nrow=10000,ncol=length(n.time)*n.sub)
  
  f.hat.final<-c()
  for(i in 1:10000){
    cor.mat<-ar.p(rho,time.pts=n.time)
    var.mat<-var.est[i]*cor.mat
    if(i < 201){resp.matrix<-rmvnorm(n.sub,mean=c(0,rep(1,length(n.time)-1)),sigma=var.mat)}
    if(i >200){resp.matrix<-rmvnorm(n.sub,mean=c(rep(0,length(n.time))),sigma=var.mat)}
    des$y<-as.vector(t(resp.matrix))
    eset[i,]<-des$y
    
    #For computational efficiency, emperical correlation function will be estimate using the first 1000 genes
    if(i<1001){
      r.raw<- residuals(lm(y~ time, data = des))
      r.frame<-data.frame(sub=des$sub,time=des$time2,resids=r.raw)
      
      
      r.wide<-acast(r.frame,sub~time,value.var="resids")
      r.cor<-rcorr(r.wide)$r
      r.tall<-data.frame(rows=rownames(r.cor)[row(r.cor)], vars=colnames(r.cor)[col(r.cor)],
                         values=c(r.cor))
      
      r.tall$d<-abs(as.numeric(as.character(r.tall$rows))-as.numeric(as.character(r.tall$vars)))
      
      f.hat<-aggregate(values~d,data=r.tall,mean)
      f.hat.final<-rbind(f.hat.final,f.hat)
      
      
    }
  }
  
  #Computing empircal correlation function and fitting the AR correlation function via Nonlinear Least Squares  
  f.bar<-aggregate(values~d,data=f.hat.final,mean)
  f.bar<-f.bar[-1,]
  f1<- values~ rho^(d)
  arfit<-nls(f1,data=f.bar,start=list(rho=f.bar$values[1]))
  
  
  csAR1<- corExp(value =
                   -1/log(coef(arfit)),
                 form = ~ time2 | sub, fixed=TRUE)
  
  csAR1 <- Initialize(csAR1, data = des)
  
  
  for (i in 1:n.genes){
    des$y<-eset[i,]
    
    
    tryCatch({
      test<-gls(model = y~ time, data = des, correlation = corExp(form = ~time2|sub),control = glsControl(opt = "optim")) 
      test2<-gls(model = y~ time, data = des, correlation = csAR1,control = glsControl(opt = "optim"))
    },error=function(e){cat("Warning: Row",i,"\n")})
    
    
    result.ar<-rbind(result.ar,coef(summary(test))[-1,4])
    result.arp<-rbind(result.arp,coef(summary(test2))[-1,4])
    
  }
  
  adjp.ar<-apply(result.ar,2,p.adjust,method="fdr")
  adjp.arp<-apply(result.arp,2,p.adjust,method="fdr")
  
  
  disc.1<-c(apply(adjp.ar,2,function(x){sum(x<.1)}),
            apply(adjp.arp,2,function(x){sum(x<.1)}))
  
  
  true.disc.1<-c(apply(adjp.ar[1:200,],2,function(x){sum(x<.1)}),
                 apply(adjp.arp[1:200,],2,function(x){sum(x<.1)}))
  
  disc<-disc+disc.1
  true.disc<-true.disc+true.disc.1
}

#Computing FDR and TPR values for plotting
fdr<-(disc-true.disc)/disc
fdrME<-1.96*sqrt(fdr*(1-fdr)/disc) 

tpr<-true.disc/n.sims/200
tprME<-1.96*sqrt(tpr*(1-tpr)/(n.sims*200))

result<-data.frame(Struc=rep(c("AR","AR Pooled"),each=4),FDR=fdr,ME=fdrME,TPR=tpr,ME_TPR=tprME)
result$time<- rep(c("T1 vs T0","T28 vs T0", "T3 vs T0", "T7 vs T0"),2)
result$time<-factor(result$time,levels=c("T1 vs T0", "T3 vs T0", "T7 vs T0","T28 vs T0"))


#Making plot

fig.6<-ggplot(data=result, aes(x=time, y=FDR, fill=Struc))+geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=FDR-ME,ymax=FDR+ME),width=.2,position=position_dodge(.9))+theme_bw()+
  scale_fill_grey(start=.4)+
  geom_hline(yintercept=0.1,size=0.75,linetype="longdash")+
  xlab("Comparison")+ylab("False Discovery Rate")+
  labs(fill='Structure') +ggtitle( bquote(AR~rho==0.2))+
  theme(legend.position="bottom",legend.background = element_rect(color = "black", 
                                                                  fill = "grey90", size = .75, linetype = "solid"),plot.title = element_text(hjust = 0.5))

fig.6

#png("Figure6.png", width = 5, height = 3.75, units="in", res=2400)
fig.6
#dev.off()








