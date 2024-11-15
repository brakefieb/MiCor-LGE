library(nlme)
library(mvtnorm)
library(ggplot2)
library(reshape2)
library(Hmisc)
library(Matrix)



#The following code directly reads in the variability estimates from the GitHub website and
#stores them. User may also download the csv file directly from the GitHub page and read in the
#file from a local directory.
var.est<-read.csv("https://raw.githubusercontent.com/brakefieb/MiCor-LGE/main/varest.csv",stringsAsFactors = TRUE)
#var.est<-read.csv("/Users/jturner/Desktop/LongitudinalProject/Codes/varest.csv",stringsAsFactors = TRUE)
var.est<-var.est[,1]


#Creating functions to generate covariance/correlation matrices for simulations.
ar.p<-function(rho,time.pts){
  t <- time.pts
  mat.final <- rho^abs(outer(t,t,"-"))
  return(mat.final)
}






#Simulating Data Set with AR (rho=0.6) structure

n.sub<-10
n.time<-c(0,1,3,7,28)

final.des<-data.frame(time=rep(n.time,n.sub),sub=paste("S",rep(1:n.sub,each=length(n.time)),sep=""))
final.des$time2<-final.des$time
final.des$time<-factor(final.des$time)

n.genes=10 
rho<-0.6


f.hat.final<-c()  #object for storing all correlations estimates for each d_ij and gene g.



set.seed(1234)
for(i in 1:10){
  cor.mat<-ar.p(rho,time.pts=n.time)
  var.mat<-var.est[i]*cor.mat
  if(i < 201){resp.matrix<-rmvnorm(n.sub,mean=c(0,rep(1,length(n.time)-1)),sigma=var.mat)}
  if(i >200){resp.matrix<-rmvnorm(n.sub,mean=c(rep(0,length(n.time))),sigma=var.mat)}
  final.des$y<-as.vector(t(resp.matrix))
  
  r.raw<- residuals(lm(y~ time, data = final.des))
  r.frame<-data.frame(sub=final.des$sub,time=final.des$time,resids=r.raw)
  
  
  r.wide<-acast(r.frame,sub~time,value.var="resids")
  r.cor<-rcorr(r.wide)$r
  
  r.tall<-data.frame(rows=rownames(r.cor)[row(r.cor)], vars=colnames(r.cor)[col(r.cor)],
                     values=c(r.cor))
  r.tall$d<-abs(as.numeric(as.character(r.tall$rows))-as.numeric(as.character(r.tall$vars)))
  
  f.hat<-aggregate(values~d,data=r.tall,function(x){tanh(mean(atanh(x)))})
  f.hat.final<-rbind(f.hat.final,f.hat)
  
}

f.hat.final<-data.frame(f.hat.final)
f.hat.final$GeneID<-factor(rep(1:10,each=11))
f.bar<-aggregate(values~d,data=f.hat.final,function(x){tanh(mean(atanh(x)))})


#Creating a new data frame of the average correlation estiamte (rho bar) along with the gene estimates (rho hat)
#For final generation of Figure 3

finalresult<-rbind(f.hat.final,data.frame(f.bar,GeneID=rep("All Genes",11)))
finalresult$Estimate<-c(rep("Ind. Genes",110),rep("All Genes",11))
finalresult$mysize<-factor(c(rep(1,110),rep(2,11)))


fig.4<-ggplot(data=finalresult,aes(x=d,y=values,group=GeneID,colour=Estimate,size=mysize))+geom_point(aes(size=Estimate))+
  geom_line(aes(linetype=Estimate,color=Estimate))+
  scale_size_manual(values=c(.5,1,1.5,.5),guide=FALSE)+scale_color_manual(values=c("#CC6666", "black"))+
  xlab(expression(d[jk]))+ylab(expression(Correlation~~f(d[jk])))

#png("Figure4.png", width = 6, height = 3.75, units="in", res=2400)
fig.4
#dev.off()



########################################################
#  Figure 5

#Removing the correlatin of 1 at d_jk=0
f.bar<-f.bar[-1,]


#Creating a starting value to estimate rho.  Using the correlation estimate at d_jk=1
init.start<-list(rho=f.bar$values[1])
arfit<-nls(values~ rho^(d),data=f.bar,start=init.start)
gausfit<-nls(values~ rho^(d^2),data=f.bar,start=init.start)


#Making predictions for plotting
predict.results<-data.frame(d=seq(0,28,.01),AR=predict(arfit,newdata=data.frame(d=seq(0,28,.01))),GAUS=predict(gausfit,newdata=data.frame(d=seq(0,28,.01))))

predict.results<-melt(predict.results,id="d")


fig.5<-  ggplot()+geom_point(data=f.bar,aes(d,values))+xlab("Lag")+ylab(expression(Correlation~~f(d[jk])))+
  geom_line(data=predict.results,aes(d,value,linetype=variable))+
  scale_linetype_manual(name="Structure",values=c("AR"=2, "GAUS"=1))


#png("Figure5.png", width = 6, height = 3.75, units="in", res=2400)
fig.5
#dev.off()






