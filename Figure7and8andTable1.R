library(nlme)
library(reshape2)
library(Hmisc)

#setwd("~/Desktop/LongitudinalProject/Codes")



#The following code directly reads in the study design matrix and the gene expression matrix 
#from the GitHub website and takes about a minute. User may also download the csv file directly 
#from the GitHub page and read in the file from a local directory to potentially save a bit of 
#time.
des<-read.csv("https://raw.githubusercontent.com/brakefieb/Misspecification-of-Correlation-Structures-in-Longitudinal-Gene-Expression-Studies-R-Scripts/main/designFluChallenge.csv",header=T,stringsAsFactors = TRUE)
#des<-read.csv("designFluChallenge.csv",header=T,stringsAsFactors = TRUE)
eset<-read.csv("https://raw.githubusercontent.com/brakefieb/Misspecification-of-Correlation-Structures-in-Longitudinal-Gene-Expression-Studies-R-Scripts/main/FluExp.csv",header=T,stringsAsFactors = TRUE)
#eset<-read.csv("FluExp.csv",header=T,stringsAsFactors = TRUE)




#it is important before running any code to check to make sure that the following are true
#correct it if it is not.
#
# 1.  To be safe make sure that eset is a matrix of numbers and that 
#  the matrix has row name id's.  Sometimes the first column will have the row ids.
# 
# 2.  After adressing one, make sure that the expression values are on the log scale.
#    Just check if there are any values in the hundreds and thousands.  If so it needs
#    to be log2 transformed.
#
# 3. Make sure the colnames of eset match up with the des$columnname in des.
# 4. Make sure that cond,time, and sub are factors in des.  Time2 should be numeric.

#Checking #1
head(eset)
#first column is the row id's so lets fix it.
rownames(eset)<-eset[,1]
eset<-as.matrix(eset[,-1])


#Checking #2
summary(as.vector(eset))
#The values are already log2 transfrmed since the max is only 15.11. If it were not
#it would be something like 2^15=32768
#If you needed to fix simply overwrite eset,  eset<-log2(eset) or log2(eset-min(eset))
#if there are negative numbers.

#Checking #3
length(setdiff(colnames(eset),as.character(des$columnname)))
#This should be 0, if not, then there are cols in eset that do not match 
#in design or vice versa.  Need to fix by subsetting the two files so that they match.

#Checking #4
is.factor(des$cond)
is.factor(des$time)
is.factor(des$sub)
is.numeric(des$time2)

#fixing variable types
table(des$time)
des$time<-factor(des$time)
des$sub<-factor(des$sub)
des$integer<-as.numeric(factor(des$time))

#Sorting desing based on subject and time point.
#First ordering everything correctly.
#This is critical to do.
index<-order(des$columnname)
des<-des[index,]
index2<-order(colnames(eset))
eset<-eset[,index2]

#Ordering by subject and time point now
index3<-order(des$sub,des$time)
des<-des[index3,]
eset<-eset[,index3]



#Initializing result to store gene level correlation estimates as a function of d_jk
f.hat.final<-c()

for(i in 1:11745){
  des$y<-eset[i,]
  r.raw<- residuals(lm(y~ time+cond+time:cond+age+gender, data = des))
  r.frame<-data.frame(sub=des$sub,time=des$time,resids=r.raw)
  
  
  r.wide<-acast(r.frame,sub~time,value.var="resids")
  r.cor<-rcorr(r.wide)$r
  r.tall<-data.frame(rows=rownames(r.cor)[row(r.cor)], vars=colnames(r.cor)[col(r.cor)],
                     values=c(r.cor))
  
  r.tall$d<-abs(as.numeric(as.character(r.tall$rows))-as.numeric(as.character(r.tall$vars)))
  
  f.hat<-aggregate(values~d,data=r.tall,mean)
  f.hat.final<-rbind(f.hat.final,f.hat)
}

#Computing emperical correlation estimates across all genes
f.bar<-aggregate(values~d,data=f.hat.final,mean)


#Plotting Figure 7  Emperical correlation function of the flu data set
f.bar<-f.bar[-1,]  #Removing d_jk=0
fig.7<-ggplot(f.bar,aes(x=d,y=values))+geom_point()+labs(x=expression(d[jk]),y=expression(rho(d[jk])))+geom_smooth()+
  ylim(0,0.75)+theme_bw()+theme(plot.title = element_text(hjust = 0.5))


#png("Figure7.png", width = 5, height = 3.75, units="in", res=2400)
fig.7
#dev.off()






#Creating Figure 8.  Fitting sphercial correlation function to the empirical estimates 
#via nonlinear least squares.


#Obtaining least squares estimates using the spherical correlation function
f1<- values~ifelse(rho1>d,(1-rho2-nugget)*(1-1.5*d/rho1+0.5*(d/rho1)^3)+rho2,rho2) 
spherfit<-nls(f1,data=f.bar,start=list(rho1=40,rho2=0.25,nugget=0.6))
coef(spherfit)


#Creating modified spherical correlation function with parameters rho1, rho2, and nugget
spher<-function(d,rho1,rho2,nugget){
  result<-c()
  for (i in 1:length(d)){
    result[i]<-ifelse(rho1>d[i],(1-rho2-nugget)*(1-1.5*d[i]/rho1+0.5*(d[i]/rho1)^3)+rho2,rho2)
  }
  return(result)
}



#png("Figure8.png", width = 5, height = 3.75, units="in", res=2400)
plot(f.bar$d,f.bar$values,ylim=c(0,.75),pch=20,ylab=expression(Correlation~~f(d[jk])),xlab=expression(d[jk]))
index<-0:120
lines(index,spher(index,rho1=36.456,rho2=0.25076,nugget=0.6513843),lwd=2)
#dev.off()





#############################################
#Table 1
#
#Applying the spherical fit to perform a genearlized linear model to determine difference over 
#time for Asymptomatic and Symptomatic groups
#
#

djk<-as.vector(dist(unique(des$time2)))
fit.vals<-spher(djk,rho1=36.456,rho2=0.25076,nugget=0.6513843)


custom.spher <- corSymm(value =fit.vals,
                        form = ~ integer | sub, fixed=TRUE)
custom.spher.init <- Initialize(custom.spher, data = des)
#corMatrix(custom.spher.init)  This line allows the user to see each subjects Correlation matrix as specified by the spherical correlation fit.




#Initializing Result objects
result.spher<-c()    #Object to store p-values for comparing HrX -Hr0 for Symp and Asymp groups using Spher
result.ar<-c()       #Object to store p-values for comparing HrX -Hr0 for Symp and Asymp groups using AR


#On rare occasions, the gls model (using REML) will fail to converge, typically when fitting AR 
#or GAUS structures. We simply note them with a warning when they occur and move forward in the
#simulation loop.
for (i in 1:11745){
  #des<-data.frame(time=rep(n.time,n.sub),sub=paste("S",rep(1:n.sub,each=length(n.time)),sep=""))
  #des$time2<-des$time
  #des$time<-factor(paste("T",des$time,sep=""))
  des$y<-eset[i,]
  #eset<-rbind(eset,des$y)
  des$cond<-factor(as.character(des$cond),levels=c("Symp","Asymp"))
  
  tryCatch({
    spher<-gls(model = y~ time+cond+time:cond+age+gender, data = des, correlation = custom.spher.init,control = glsControl(opt = "optim")) 
    ar<-gls(model = y~ time+cond+time:cond+age+gender, data = des, correlation = corExp(form = ~time2|sub),control = glsControl(opt = "optim"))
  },error=function(e){cat("Warning: Row",i,"\n")})
  
  des$cond<-factor(as.character(des$cond),levels=c("Asymp","Symp"))
  tryCatch({
    spher2<-gls(model = y~ time+cond+time:cond+age+gender, data = des, correlation = custom.spher.init,control = glsControl(opt = "optim")) 
    ar2<-gls(model = y~ time+cond+time:cond+age+gender, data = des, correlation = corExp(form = ~time2|sub),control = glsControl(opt = "optim"))
  },error=function(e){cat("Warning: Row",i,"\n")})
  
  
  #Combining the HrX vs Hr0 p-value comparisons for Symptomatic (spher and ar) with the same comparisons for Asymptomatic (spher2 and ar2)
  result.spher<-rbind(result.spher,c(coef(summary(spher))[3:16,4],coef(summary(spher2))[3:16,4]))
  result.ar<-rbind(result.ar,c(coef(summary(ar))[3:16,4],coef(summary(ar2))[3:16,4]))
  
}

adjp.spher<-apply(result.spher,2,p.adjust,method="fdr")
adjp.ar<-apply(result.ar,2,p.adjust,method="fdr")

rejected.spher<-apply(adjp.spher<.1,2,sum)
rejected.ar<-apply(adjp.ar<.1,2,sum)

table.1<-data.frame(Comparison=paste("HR",unique(des$time2)[-(1:2)],"-Base",sep=""),
                    SympSPHER=rejected.spher[1:14],
                    SympAR=rejected.ar[1:14],
                    AsympSPH=rejected.spher[15:28],
                    AsympAR=rejected.ar[15:28])

View(table.1)

write.csv(table.1,file="Table1.csv")


#For curious readers and reviewers.  The significant genes in the Spher model are almost always contined within
#the gene lists of the AR model.  The Spher model just has additional significant genes.
counts.spher<- adjp.spher<.1
counts.ar<- adjp.ar<.1

common<-c()
for (i in 1:26){
  common[i]<-length(intersect(which(counts.spher[,i]==TRUE), which(counts.ar[,i]==TRUE)))
}

CommonGenes<-data.frame(Comparison=paste("HR",unique(des$time2)[-(1:2)],"-Base",sep=""),
                        Symp=common[1:14],
                        Asymp=common[15:28])

View(CommonGenes)


write.csv(CommonGenes,file="CommonGenes.csv")












