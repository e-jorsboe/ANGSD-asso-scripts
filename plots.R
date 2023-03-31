
cex<-1.5

###########################################################

param <- list(sd=1,Nind=1000,beta=0.6,freq=c(0.9,0.1),delta=3,simulations=10000)
param$Q=rbind(rep(c(0.05,0.95,0.5,0.95),each=param$Nind/4),1-rep(c(0.05,0.95,0.5,0.95),each=param$Nind/4))
param$fInd <- colSums(param$Q*param$freq)
param$depthGroup<- c(4,1)
param$threshold <- 1e-3

deltaValue<-5
gammaValue<-1


## remember power is rm 0 reads, power2 is keep 0 reads
load("Fig2.Rdata")

###########################################################

pdf("Fig2a.pdf")
coll <- c("darkred","darkblue","darkgreen","red","blue")

ymax<-max(c(sapply(power2,function(x) x$powF)/param$threshold,sapply(power2,function(x) x$powInd)/param$threshold,sapply(power,function(x) x$powG)/param$threshold,
            sapply(power2,function(x) x$powDosF)/param$threshold,sapply(power2,function(x) x$powDosInd)/param$threshold))

## round up to nearst 10
ymax<-ceiling(ymax/10)*10

plot(range,sapply(power2,function(x) x$powF)/param$threshold,type="l",col=coll[1],lwd=3,ylab="false pos. rate / expected false pos. rate",xlab="delta",
     main="false postive rate enrichment",cex.axis = cex,cex.lab = cex,ylim=c(0,50))
lines(range,sapply(power2,function(x) x$powInd)/param$threshold,col=coll[2],lwd=3)
lines(range,sapply(power2,function(x) x$powG)/param$threshold,col=coll[3],lwd=3)
lines(range,sapply(power2,function(x) x$powDosF)/param$threshold,col=coll[4],lwd=3,lty=2)
lines(range,sapply(power2,function(x) x$powDosInd)/param$threshold,col=coll[5],lwd=3,lty=2)


legend("topleft",col=coll,c("ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype",expression(paste("dosages (", f,")",sep="")),expression(paste("dosages (", pi,")",sep=""))),cex=cex,pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,2,2),lwd=2,bty = "n")


dev.off()



cex<-1.4

pdf("Fig2b.pdf")
coll <- c("darkred","darkblue")
plot(range,sapply(power2,function(x) x$diffF),type="l",col=coll[1],lwd=3,ylab="bias",xlab="delta",main="bias of effect size",
     ylim=c(min(c(sapply(power2,function(x) x$diffF),sapply(power2,function(x) x$diffInd))),max(c(sapply(power2,function(x) x$diffF),sapply(power2,function(x) x$diffInd)))),cex.axis = cex,cex.lab = cex)
lines(range,sapply(power2,function(x) x$diffInd),col=coll[2],lwd=3)
legend("bottomleft",col=coll,c("bias ANGSD-asso (f)",expression(paste("bias ANGSD-asso (", pi,")",sep=""))),cex = cex,bty = "n",pch=c(NA,NA),lty=c(1,1),lwd=2)
dev.off()

cex<-1.5



pdf("Fig2c.pdf",width=20)
cex<-2

par(mfrow=c(1,2))
par(mar=c(5.1,5.1,4.1,2.1))
b<-barplot(param$Q[1:2,],ylab="admix prop. population 1",xlab="individuals",cex.axis = cex,cex.lab = cex,ylim=c(0,1),main="admixture proportions",
           cex.main=cex,space=0,border=NA,col=c("black","grey"))

abline(v=0)
abline(v=250)
abline(v=500)
abline(v=750)
abline(v=1000)
segments(0, 1, 1000, 1)

mtext(line=1,paste0("N=",param$Nind/4),1,at = mean(c(quantile(b)[1],quantile(b)[2])),cex=cex)
mtext(line=1,paste0("N=",param$Nind/4),1,at = mean(c(quantile(b)[2],quantile(b)[3])),cex=cex)
mtext(line=1,paste0("N=",param$Nind/4),1,at = mean(c(quantile(b)[3],quantile(b)[4])),cex=cex)
mtext(line=1,paste0("N=",param$Nind/4),1,at = mean(c(quantile(b)[4],quantile(b)[5])),cex=cex)

b<-barplot(param$depthGroup,ylab="avg. depth",xlab="individuals",cex.axis = cex,cex.lab = cex,main="sequencing depth",cex.main=cex)
mtext(line=1,paste0("N=",param$Nind/2),1,at = mean(c(quantile(b)[1],quantile(b)[3])),cex=cex)
mtext(line=1,paste0("N=",param$Nind/2),1,at = mean(c(quantile(b)[3],quantile(b)[5])),cex=cex)

dev.off()

cex<-1.5



######################################
## real power
#########################################


## remember power is rm 0 reads, power2 is keep 0 reads
load("Fig3.Rdata")


cex<-1.4


pdf("Fig3a.pdf")
coll <- c("darkred","darkblue","darkgreen","red","blue")

plot(range,sapply(power2,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="beta",
     main="",
     ylim=c(0,1),cex.axis = cex,cex.lab = cex)
lines(range,sapply(power2,function(x) x$powInd),col=coll[2],lwd=3)
lines(range,sapply(power2,function(x) x$powG),col=coll[3],lwd=3)
lines(range,sapply(power2,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
lines(range,sapply(power2,function(x) x$powDosInd),col=coll[5],lwd=3,lty=2)


legend("topleft",col=coll,c("ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype",expression(paste("dosages (", f,")",sep="")),expression(paste("dosages (", pi,")",sep=""))),cex=(cex-0.1),pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,2,2),lwd=2,bty = "n")


dev.off()

pdf("Fig3b.pdf")
coll <- c("darkred","darkblue")
plot(range,sapply(power2,function(x) x$diffF),type="l",col=coll[1],lwd=3,ylab="bias",xlab="beta",main="bias of effect size",
     ylim=c(min(c(sapply(power2,function(x) x$diffF),sapply(power2,function(x) x$diffInd))),max(c(sapply(power2,function(x) x$diffF),sapply(power2,function(x) x$diffInd)))),cex.axis = cex,cex.lab = cex)
lines(range,sapply(power2,function(x) x$diffInd),col=coll[2],lwd=3)
legend("bottomleft",col=coll,c("bias ANGSD-asso (f)",expression(paste("bias ANGSD-asso (", pi,")",sep=""))),cex = cex,bty = "n",pch=c(NA,NA),lty=c(1,1),lwd=2)



dev.off()



########################################
##### SCENARIOS WITHOUT ADMIXTURE
########################################


param <- list(sd=1,Nind=1000,beta=0.6,freq=c(0.45,0.45),delta=3,simulations=10000)
param$Q=rbind(rep(c(1,1,1,1),each=param$Nind/4),(1-rep(c(1,1,1,1),each=param$Nind/4)))
param$fInd <- colSums(param$Q*param$freq)
param$depthGroup<- c(4,1)
param$threshold <- 1e-3

deltaValue<-5
gammaValue<-0


## remember power is rm 0 reads, power2 is keep 0 reads
load("SuppFig1.Rdata")


###########################
## effect of depth
#################################


pdf("SuppFig1a.pdf")
coll <- c("darkred","darkblue","darkgreen","red","blue")


plot(range,sapply(power2,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="false positive rate",xlab="delta",
     main="",
     ylim=c(0,0.2),cex.axis = cex,cex.lab = cex)
lines(range,sapply(power2,function(x) x$powInd),col=coll[2],lwd=3)
lines(range,sapply(power2,function(x) x$powG),col=coll[3],lwd=3)
lines(range,sapply(power2,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
lines(range,sapply(power2,function(x) x$powDosInd),col=coll[5],lwd=3,lty=2)


legend("topleft",col=coll,c("ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype",expression(paste("dosages (", f,")",sep="")),expression(paste("dosages (", pi,")",sep=""))),cex=cex,pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,2,2),lwd=2,bty = "n")


dev.off()


pdf("SuppFig1b.pdf")
coll <- c("darkred","darkblue")
plot(range,sapply(power2,function(x) x$diffF),type="l",col=coll[1],lwd=3,ylab="bias",xlab="delta",main="bias of effect size",cex.axis = cex,cex.lab = cex,ylim=c(-0.05,0.05))
lines(range,sapply(power2,function(x) x$diffInd),col=coll[2],lwd=3)
legend("topleft",col=coll,c("bias ANGSD-asso (f)",expression(paste("bias ANGSD-asso (", pi,")",sep=""))),cex = cex,bty = "n",pch=c(NA,NA),lty=c(1,1),lwd=2)
dev.off()

pdf("SuppFig1c.pdf")
coll <- c("darkred","darkblue","darkgreen","red","blue")

plot(range,sapply(power2,function(x) x$powF)/param$threshold,type="l",col=coll[1],lwd=3,ylab="false pos. rate / expected false pos. rate",xlab="delta",
     main="false postive rate enrichment",cex.axis = cex,cex.lab = cex,ylim=c(0,50))
lines(range,sapply(power2,function(x) x$powInd)/param$threshold,col=coll[2],lwd=3)
lines(range,sapply(power2,function(x) x$powG)/param$threshold,col=coll[3],lwd=3)
lines(range,sapply(power2,function(x) x$powDosF)/param$threshold,col=coll[4],lwd=3,lty=2)
lines(range,sapply(power2,function(x) x$powDosInd)/param$threshold,col=coll[5],lwd=3,lty=2)


legend("topleft",col=coll,c("ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype",expression(paste("dosages (", f,")",sep="")),expression(paste("dosages (", pi,")",sep=""))),cex=cex,pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,2,2),lwd=2,bty = "n")


dev.off()



cex<-1.5

pdf("SuppFig1c.pdf",width=20)
cex<-2
par(mfrow=c(1,2))
par(mar=c(5.1,5.1,4.1,2.1))
b<-barplot(param$Q[1:2,],ylab="admix prop. population 1",xlab="individuals",cex.axis = cex,cex.lab = cex,
           ylim=c(0,1),main="admixture proportions",cex.main=cex,space=0,border=NA,col=c("black","grey"))


abline(v=0)
abline(v=250)
abline(v=500)
abline(v=750)
abline(v=1000)
segments(0, 1, 1000, 1)

mtext(line=1,paste0("N=",table(param$Q[1,])[1]),1,at = mean(c(quantile(b)[1],quantile(b)[5])),cex=cex)


b<-barplot(param$depthGroup,ylab="avg. depth",xlab="individuals",cex.axis = cex,cex.lab = cex,main="sequencing depth",cex.main=cex)


mtext(line=1,paste0("N=",param$Nind/2),1,at = mean(c(quantile(b)[1],quantile(b)[3])),cex=cex)
mtext(line=1,paste0("N=",param$Nind/2),1,at = mean(c(quantile(b)[3],quantile(b)[5])),cex=cex)


dev.off()

cex<-1.5


###################################################################

## remember power is rm 0 reads, power2 is keep 0 reads
load("SuppFig2.Rdata")

cex<-1.4


pdf("SuppFig2a.pdf")
coll <- c("darkred","darkblue","darkgreen","red","blue")

plot(range,sapply(power2,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="beta",
     main="statistical power",cex.axis = cex,cex.lab = cex,ylim=c(0,1))
lines(range,sapply(power2,function(x) x$powInd),col=coll[2],lwd=3)
lines(range,sapply(power2,function(x) x$powG),col=coll[3],lwd=3)
lines(range,sapply(power2,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
lines(range,sapply(power2,function(x) x$powDosInd),col=coll[5],lwd=3,lty=2)

legend("topleft",col=coll,c("ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype",expression(paste("dosages (", f,")",sep="")),expression(paste("dosages (", pi,")",sep=""))),cex=(cex-0.1),pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,2,2),lwd=2,bty = "n")

dev.off()

pdf("SuppFig2b.pdf")
coll <- c("darkred","darkblue")

plot(range,sapply(power2,function(x) x$diffF),type="l",col=coll[1],lwd=3,ylab="bias",xlab="beta",main="bias of effect size",cex.axis = cex,cex.lab = cex,ylim=c(-0.05,0.05))
lines(range,sapply(power2,function(x) x$diffInd),col=coll[2],lwd=3)
legend("topleft",col=coll,c("bias ANGSD-asso (f)",expression(paste("bias ANGSD-asso (", pi,")",sep=""))),cex = cex,bty = "n",pch=c(NA,NA),lty=c(1,1),lwd=2)

dev.off()



############################################################
############################################################

cex<-1.5


load("SuppFig3and4.Rdata")

coll <- c("darkred","darkblue","darkgreen","red","blue")
ylim<-c(0,1)


pdf("SuppFig3a.pdf")

plot(range,sapply(power,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="effect size",
     main="remove 0 reads, est. freq.",cex.axis = cex,cex.lab = cex,ylim=ylim)
lines(range,sapply(power,function(x) x$powG),col=coll[3],lwd=3)
lines(range,sapply(power,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
legend("topleft",col=coll[c(1,3,4)],c("ANGSD-asso (f)","true genotype","dosage (f)"),cex=(cex-0.1),pch=c(NA,NA,NA),lty=c(1,1,2),lwd=2,bty = "n")

dev.off()

pdf("SuppFig3b.pdf")

plot(range,sapply(power2,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="effect size",
     main="keep 0 reads, est. freq.",cex.axis = cex,cex.lab = cex,ylim=ylim)
lines(range,sapply(power2,function(x) x$powG),col=coll[3],lwd=3)
lines(range,sapply(power2,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
legend("topleft",col=coll[c(1,3,4)],c("ANGSD-asso (f)","true genotype","dosage (f)"),cex=(cex-0.1),pch=c(NA,NA,NA),lty=c(1,1,2),lwd=2,bty = "n")

dev.off()

pdf("SuppFig3c.pdf")

plot(range,sapply(power3,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="effect size",
     main="remove 0 reads, oracle freq.",cex.axis = cex,cex.lab = cex,ylim=ylim)
lines(range,sapply(power3,function(x) x$powG),col=coll[3],lwd=3)
lines(range,sapply(power3,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
legend("topleft",col=coll[c(1,3,4)],c("ANGSD-asso (f)","true genotype","dosage (f)"),cex=(cex-0.1),pch=c(NA,NA,NA),lty=c(1,1,2),lwd=2,bty = "n")

dev.off()

pdf("SuppFig3d.pdf")

plot(range,sapply(power4,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="effect size",
     main="keep 0 reads, oracle freq.",cex.axis = cex,cex.lab = cex,ylim=ylim)
lines(range,sapply(power4,function(x) x$powG),col=coll[3],lwd=3)
lines(range,sapply(power4,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
legend("topleft",col=coll[c(1,3,4)],c("ANGSD-asso (f)","true genotype","dosage (f)"),cex=(cex-0.1),pch=c(NA,NA,NA),lty=c(1,1,2),lwd=2,bty = "n")

dev.off()



coll <- c("darkred","darkblue","darkgreen","red","blue")
ylim<-c(-0.2,0.2)


pdf("SuppFig4a.pdf")

plot(range,sapply(power,function(x) x$diffF),type="l",col=coll[1],lwd=3,ylab="bias",xlab="effect size",main="remove 0 reads, est. freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(range,sapply(power,function(x) x$diffDosF),col=coll[4],lwd=3,lty=2)
lines(range,sapply(power,function(x) x$diffGeno),col=coll[3],lwd=3)
abline(h=0)
legend("topright",col=c(coll[1],coll[4],coll[3]),c("bias ANGSD-asso (f)","bias dosage (f)","bias genotype"),cex = cex,bty = "n",pch=c(NA,NA,NA),lty=c(1,2,1),lwd=2)

dev.off()

pdf("SuppFig4b.pdf")

plot(range,sapply(power2,function(x) x$diffF),type="l",col=coll[1],lwd=3,ylab="bias",xlab="effect size",main="keep 0 reads, est. freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(range,sapply(power2,function(x) x$diffDosF),col=coll[4],lwd=3,lty=2)
lines(range,sapply(power2,function(x) x$diffGeno),col=coll[3],lwd=3)
abline(h=0)
legend("topright",col=c(coll[1],coll[4],coll[3]),c("bias ANGSD-asso (f)","bias dosage (f)","bias genotype"),cex = cex,bty = "n",pch=c(NA,NA,NA),lty=c(1,2,1),lwd=2)

dev.off()

pdf("SuppFig4c.pdf")

plot(range,sapply(power3,function(x) x$diffF),type="l",col=coll[1],lwd=3,ylab="bias",xlab="effect size",main="remove 0 reads, oracle freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(range,sapply(power3,function(x) x$diffDosF),col=coll[4],lwd=3,lty=2)
lines(range,sapply(power3,function(x) x$diffGeno),col=coll[3],lwd=3)
abline(h=0)
legend("topright",col=c(coll[1],coll[4],coll[3]),c("bias ANGSD-asso (f)","bias dosage (f)","bias genotype"),cex = cex,bty = "n",pch=c(NA,NA,NA),lty=c(1,2,1),lwd=2)

dev.off()

pdf("SuppFig4d.pdf")

plot(range,sapply(power4,function(x) x$diffF),type="l",col=coll[1],lwd=3,ylab="bias",xlab="effect size",main="keep 0 reads, oracle freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(range,sapply(power4,function(x) x$diffDosF),col=coll[4],lwd=3,lty=2)
lines(range,sapply(power4,function(x) x$diffGeno),col=coll[3],lwd=3)
abline(h=0)
legend("topright",col=c(coll[1],coll[4],coll[3]),c("bias ANGSD-asso (f)","bias dosage (f)","bias genotype"),cex = cex,bty = "n",pch=c(NA,NA,NA),lty=c(1,2,1),lwd=2)

dev.off()


###############################################
###############################################

load("SuppFig5and6.Rdata")

coll <- c("darkred","darkblue","darkgreen","red","blue")
ylim=c(0,1)


pdf("SuppFig5a.pdf")

plot(range,sapply(power,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="beta",
     main="remove 0 reads, est. freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(range,sapply(power,function(x) x$powInd),col=coll[2],lwd=3)
lines(range,sapply(power,function(x) x$powG),col=coll[3],lwd=3)
lines(range,sapply(power,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
lines(range,sapply(power,function(x) x$powDosInd),col=coll[5],lwd=3,lty=2)
legend("topleft",col=coll,c("ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype",expression(paste("dosage (", f,")",sep="")),expression(paste("dosage (", pi,")",sep=""))),cex=(cex-0.1),pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,2,2),lwd=2,bty = "n")

dev.off()

pdf("SuppFig5b.pdf")

plot(range,sapply(power2,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="beta",
     main="keep 0 reads, est. freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(range,sapply(power2,function(x) x$powInd),col=coll[2],lwd=3)
lines(range,sapply(power2,function(x) x$powG),col=coll[3],lwd=3)
lines(range,sapply(power2,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
lines(range,sapply(power2,function(x) x$powDosInd),col=coll[5],lwd=3,lty=2)
legend("topleft",col=coll,c("ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype",expression(paste("dosage (", f,")",sep="")),expression(paste("dosage (", pi,")",sep=""))),cex=(cex-0.1),pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,2,2),lwd=2,bty = "n")

dev.off()

pdf("SuppFig5c.pdf")

plot(range,sapply(power3,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="beta",
     main="remove 0 reads, oracle freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(range,sapply(power3,function(x) x$powInd),col=coll[2],lwd=3)
lines(range,sapply(power3,function(x) x$powG),col=coll[3],lwd=3)
lines(range,sapply(power3,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
lines(range,sapply(power3,function(x) x$powDosInd),col=coll[5],lwd=3,lty=2)
legend("topleft",col=coll,c("ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype",expression(paste("dosage (", f,")",sep="")),expression(paste("dosage (", pi,")",sep=""))),cex=(cex-0.1),pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,2,2),lwd=2,bty = "n")

dev.off()

pdf("SuppFig5d.pdf")

plot(range,sapply(power4,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="beta",
     main="keep 0 reads, oracle freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(range,sapply(power4,function(x) x$powInd),col=coll[2],lwd=3)
lines(range,sapply(power4,function(x) x$powG),col=coll[3],lwd=3)
lines(range,sapply(power4,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
lines(range,sapply(power4,function(x) x$powDosInd),col=coll[5],lwd=3,lty=2)
legend("topleft",col=coll,c("ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype",expression(paste("dosage (", f,")",sep="")),expression(paste("dosage (", pi,")",sep=""))),cex=(cex-0.1),pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,2,2),lwd=2,bty = "n")

dev.off()

coll <- c("darkblue","darkblue","darkgreen")
ylim=c(-0.2,0.2)

pdf("SuppFig6a.pdf")

plot(range,sapply(power,function(x) x$diffDosInd),type="l",col=coll[1],lty=2,lwd=3,ylab="bias",xlab="effect size",main="remove 0 reads, est. freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(range,sapply(power,function(x) x$diffInd),col=coll[2],lwd=3)
lines(range,sapply(power,function(x) x$diffGeno),col=coll[3],lwd=3)
abline(h=0)
legend("top",col=coll,c(expression(paste("ANGSD-asso (", pi,")",sep="")),expression(paste("bias dosage (", pi,")",sep="")),"bias genotype"),cex = cex,bty = "n",pch=c(NA,NA,NA),lty=c(1,2,1),lwd=2)

dev.off()

pdf("SuppFig6b.pdf")

plot(range,sapply(power2,function(x) x$diffDosInd),type="l",col=coll[1],lty=2,lwd=3,ylab="bias",xlab="effect size",main="keep 0 reads, est. freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(range,sapply(power2,function(x) x$diffInd),col=coll[2],lwd=3)
lines(range,sapply(power2,function(x) x$diffGeno),col=coll[3],lwd=3)
abline(h=0)
legend("top",col=coll,c(expression(paste("ANGSD-asso (", pi,")",sep="")),expression(paste("bias dosage (", pi,")",sep="")),"bias genotype"),cex = cex,bty = "n",pch=c(NA,NA,NA),lty=c(1,2,1),lwd=2)

dev.off()

pdf("SuppFig6c.pdf")

plot(range,sapply(power3,function(x) x$diffDosInd),type="l",col=coll[1],lty=2,lwd=3,ylab="bias",xlab="effect size",main="remove 0 reads, oracle freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(range,sapply(power3,function(x) x$diffInd),col=coll[2],lwd=3)
lines(range,sapply(power3,function(x) x$diffGeno),col=coll[3],lwd=3)
abline(h=0)
legend("top",col=coll,c(expression(paste("ANGSD-asso (", pi,")",sep="")),expression(paste("bias dosage (", pi,")",sep="")),"bias genotype"),cex = cex,bty = "n",pch=c(NA,NA,NA),lty=c(1,2,1),lwd=2)

dev.off()

pdf("SuppFig6d.pdf")

plot(range,sapply(power4,function(x) x$diffDosInd),type="l",col=coll[1],lty=2,lwd=3,ylab="bias",xlab="effect size",main="keep 0 reads, oracle freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(range,sapply(power4,function(x) x$diffInd),col=coll[2],lwd=3)
lines(range,sapply(power4,function(x) x$diffGeno),col=coll[3],lwd=3)
abline(h=0)
legend("top",col=coll,c(expression(paste("ANGSD-asso (", pi,")",sep="")),expression(paste("bias dosage (", pi,")",sep="")),"bias genotype"),cex = cex,bty = "n",pch=c(NA,NA,NA),lty=c(1,2,1),lwd=2)

dev.off()


##########################################################
##########################################################



load("SuppFig7and8.Rdata")

coll <- c("darkred","darkblue","darkgreen","red","blue")
ylim=c(0,1)

pdf("SuppFig7a.pdf")

plot(exp(range),sapply(power,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="OR",
     main="remove 0 reads, est. freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(exp(range),sapply(power,function(x) x$powInd),col=coll[2],lwd=3)
lines(exp(range),sapply(power,function(x) x$powG),col=coll[3],lwd=3)
lines(exp(range),sapply(power,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
lines(exp(range),sapply(power,function(x) x$powDosInd),col=coll[5],lwd=3,lty=2)
legend("topleft",col=coll,c("ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype",expression(paste("dosage (", f,")",sep="")),expression(paste("dosage (", pi,")",sep=""))),cex=(cex-0.1),pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,2,2),lwd=2,bty = "n")

dev.off()

pdf("SuppFig7b.pdf")

plot(exp(range),sapply(power2,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="OR",
     main="keep 0 reads, est. freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(exp(range),sapply(power2,function(x) x$powInd),col=coll[2],lwd=3)
lines(exp(range),sapply(power2,function(x) x$powG),col=coll[3],lwd=3)
lines(exp(range),sapply(power2,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
lines(exp(range),sapply(power2,function(x) x$powDosInd),col=coll[5],lwd=3,lty=2)
legend("topleft",col=coll,c("ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype",expression(paste("dosage (", f,")",sep="")),expression(paste("dosage (", pi,")",sep=""))),cex=(cex-0.1),pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,2,2),lwd=2,bty = "n")

dev.off()

pdf("SuppFig7c.pdf")

plot(exp(range),sapply(power3,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="OR",
     main="remove 0 reads, oracle freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(exp(range),sapply(power3,function(x) x$powInd),col=coll[2],lwd=3)
lines(exp(range),sapply(power3,function(x) x$powG),col=coll[3],lwd=3)
lines(exp(range),sapply(power3,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
lines(exp(range),sapply(power3,function(x) x$powDosInd),col=coll[5],lwd=3,lty=2)
legend("topleft",col=coll,c("ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype",expression(paste("dosage (", f,")",sep="")),expression(paste("dosage (", pi,")",sep=""))),cex=(cex-0.1),pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,2,2),lwd=2,bty = "n")

dev.off()

pdf("SuppFig7d.pdf")

plot(exp(range),sapply(power4,function(x) x$powF),type="l",col=coll[1],lwd=3,ylab="power",xlab="OR",
     main="keep 0 reads, oracle freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(exp(range),sapply(power4,function(x) x$powInd),col=coll[2],lwd=3)
lines(exp(range),sapply(power4,function(x) x$powG),col=coll[3],lwd=3)
lines(exp(range),sapply(power4,function(x) x$powDosF),col=coll[4],lwd=3,lty=2)
lines(exp(range),sapply(power4,function(x) x$powDosInd),col=coll[5],lwd=3,lty=2)
legend("topleft",col=coll,c("ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype",expression(paste("dosage (", f,")",sep="")),expression(paste("dosage (", pi,")",sep=""))),cex=(cex-0.1),pch=c(NA,NA,NA,NA,NA),lty=c(1,1,1,2,2),lwd=2,bty = "n")

dev.off()

coll <- c("darkblue","darkblue","darkgreen")
ylim=c(-0.2,0.2)

pdf("SuppFig8a.pdf")

plot(exp(range),sapply(power,function(x) x$diffDosInd),type="l",col=coll[2],lty=2,lwd=3,ylab="bias beta",xlab="OR",main="remove 0 reads, est. freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(exp(range),sapply(power,function(x) x$diffInd),col=coll[2],lwd=3)
lines(exp(range),sapply(power,function(x) x$diffGeno),col=coll[3],lwd=3)
abline(h=0)
legend("top",col=coll,c(expression(paste("bias dosage (", pi,")",sep="")),expression(paste("bias ANGSD-asso (", pi,")",sep="")),"bias genotype"),cex = cex,bty = "n",pch=c(NA,NA,NA),lty=c(2,1,1),lwd=2)

dev.off()

pdf("SuppFig8b.pdf")

plot(exp(range),sapply(power2,function(x) x$diffDosInd),type="l",col=coll[2],lty=2,lwd=3,ylab="bias beta",xlab="OR",main="keep 0 reads, est. freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(exp(range),sapply(power2,function(x) x$diffInd),col=coll[2],lwd=3)
lines(exp(range),sapply(power2,function(x) x$diffGeno),col=coll[3],lwd=3)
abline(h=0)
legend("top",col=coll,c(expression(paste("bias dosage (", pi,")",sep="")),expression(paste("bias ANGSD-asso (", pi,")",sep="")),"bias genotype"),cex = cex,bty = "n",pch=c(NA,NA,NA),lty=c(2,1,1),lwd=2)

dev.off()

pdf("SuppFig8c.pdf")

plot(exp(range),sapply(power3,function(x) x$diffDosInd),type="l",col=coll[2],lty=2,lwd=3,ylab="bias beta",xlab="OR",main="remove 0 reads, oracle freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(exp(range),sapply(power3,function(x) x$diffInd),col=coll[2],lwd=3)
lines(exp(range),sapply(power3,function(x) x$diffGeno),col=coll[3],lwd=3)
abline(h=0)
legend("top",col=coll,c(expression(paste("bias dosage (", pi,")",sep="")),expression(paste("bias ANGSD-asso (", pi,")",sep="")),"bias genotype"),cex = cex,bty = "n",pch=c(NA,NA,NA),lty=c(2,1,1),lwd=2)

dev.off()

pdf("SuppFig8d.pdf")

plot(exp(range),sapply(power4,function(x) x$diffDosInd),type="l",col=coll[2],lty=2,lwd=3,ylab="bias beta",xlab="OR",main="keep 0 reads, oracle freq.",
     ylim=ylim,cex.axis = cex,cex.lab = cex)
lines(exp(range),sapply(power4,function(x) x$diffInd),col=coll[2],lwd=3)
lines(exp(range),sapply(power4,function(x) x$diffGeno),col=coll[3],lwd=3)
abline(h=0)
legend("top",col=coll,c(expression(paste("bias dosage (", pi,")",sep="")),expression(paste("bias ANGSD-asso (", pi,")",sep="")),"bias genotype"),cex = cex,bty = "n",pch=c(NA,NA,NA),lty=c(2,1,1),lwd=2)

dev.off()


############################################################
############################################################


biasLim<-c(-0.4,0.4)


load("SuppFig12.Rdata")


pdf("SuppFig12a.pdf")
plot(betaRange[1:13],snptestVectorP[1:13],type="l",ylim=c(0,1),ylab="power",xlab="beta",col="black",lwd=3,cex.axis = cex,cex.lab = cex)
points(betaRange[1:13],snptestVectorP2[1:13],type="l",lwd=3,col="black",lty=2)
points(betaRange[1:13],angsdVectorP[1:13],type="l",lwd=3,col="darkgrey")
points(betaRange[1:13],angsdVectorP2[1:13],type="l",lwd=3,col="darkgrey",lty=2)
points(betaRange[1:13],genoVectorP[1:13],type="l",col="darkgreen",lwd=3)
cols<-c("black","black","darkgrey","darkgrey","darkgreen")
legend("topleft",c("SNPTEST (f)",expression(paste("SNPTEST (", pi,")",sep="")),"ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype"),col=cols,cex=(cex-0.1),pch=c(NA,NA,NA,NA,NA),lty=c(1,2,1,2,1),lwd=2,bty = "n")
dev.off()

pdf("SuppFig12b.pdf")
plot(betaRange[1:13],snptestVector[1:13],type="l",ylim=biasLim,ylab="bias",xlab="beta",col="black",lwd=3,cex.axis = cex,cex.lab = cex)
points(betaRange[1:13],genoVector[1:13],type="l",lwd=3,col="darkgreen")
points(betaRange[1:13],snptestVector2[1:13],type="l",lwd=3,col="black",lty=2)
points(betaRange[1:13],angsdVector[1:13],type="l",lwd=3,col="darkgrey")
points(betaRange[1:13],angsdVector2[1:13],type="l",lwd=3,col="darkgrey",lty=2)
cols<-c("black","black","darkgrey","darkgrey","darkgreen")
legend("topleft",c("SNPTEST (f)",expression(paste("SNPTEST (", pi,")",sep="")),"ANGSD-asso (f)",expression(paste("ANGSD-asso (", pi,")",sep="")),"true genotype"),col=cols,cex=(cex-0.1),pch=c(NA,NA,NA,NA,NA),lty=c(1,2,1,2,1),lwd=2,bty = "n")
dev.off()
