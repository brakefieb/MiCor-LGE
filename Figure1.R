library(ggplot2)
library(reshape2)

#setwd("~/Desktop/LongitudinalProject/Codes")

#Figure 1
d<-seq(0,5,.01)
rho<-0.6
cors<-data.frame(d=d,CS=rep(rho,length(d)),AR=rho^d,GAUS=rho^(d^2))

cors.long<-melt(cors,id="d")
names(cors.long)<-c("d","Structure","Correlation")

#png("Figure1.png", width = 5, height = 3.75, units="in", res=4800)
ggplot(data=cors.long, aes(x=d, y=Correlation, group=Structure)) +theme_bw()+
  geom_line(aes(linetype=Structure),size=1.25)+
  scale_linetype_manual(values=c("solid","twodash", "dotted"))+
  xlab(expression(d[jk]))+ylab(expression(f(d[jk])))
#dev.off()
