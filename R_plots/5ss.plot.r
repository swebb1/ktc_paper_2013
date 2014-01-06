##Read in tables
p0<-read.table("Galaxy216-[p0.5ss.windows].interval")
p30<-read.table("Galaxy228-[p30.5ss.windows].interval")
wt0<-read.table("Galaxy240-[wt0.5ss.windows].interval")
wt30<-read.table("Galaxy252-[wt30.5ss.windows].interval")
p0_2<-read.table("Galaxy264-[p0_200.5ss.windows].interval")
p30_2<-read.table("Galaxy276-[p30_200.5ss.windows].interval")
wt0_2<-read.table("Galaxy288-[wt0_200.5ss.windows].interval")
wt30_2<-read.table("Galaxy300-[wt30_200.5ss.windows].interval")

#Get table dimensions
l<-dim(p0)[2]
g<-dim(p0)[1]

#Remove 0 values
sum<-vector()
for(i in 1:g){sum<-append(sum,sum(p0[i,4:l])+sum(p30[i,4:l])+sum(wt0[i,4:l])+sum(wt30[i,4:l]))}
p0<-p0[which(sum>1),]
p30<-p30[which(sum>1),]
wt0<-wt0[which(sum>1),]
wt30<-wt30[which(sum>1),]

g1<-dim(p0)[1]

sum<-vector()
for(i in 1:g){sum<-append(sum,sum(p0_2[i,4:l])+sum(p30_2[i,4:l])+sum(wt0_2[i,4:l])+sum(wt30_2[i,4:l]))}
p0_2<-p0_2[which(sum>1),]
p30_2<-p30_2[which(sum>1),]
wt0_2<-wt0_2[which(sum>1),]
wt30_2<-wt30_2[which(sum>1),]

g2<-dim(p0_2)[1]

#write.table(p0[,1:3],file="247_5ss_genes",sep="\t")
#write.table(p0_2[,1:3],file="249_5ss_genes",sep="\t")

#Get mean values
tp0<-vector()
tp30<-vector()
twt0<-vector()
twt30<-vector()
tp0_2<-vector()
tp30_2<-vector()
twt0_2<-vector()
twt30_2<-vector()


for(i in 4:l){tp0<-append(tp0,mean(p0[,i]))}
for(i in 4:l){tp30<-append(tp30,mean(p30[,i]))}
for(i in 4:l){twt0<-append(twt0,mean(wt0[,i]))}
for(i in 4:l){twt30<-append(twt30,mean(wt30[,i]))}
for(i in 4:l){tp0_2<-append(tp0_2,mean(p0_2[,i]))}
for(i in 4:l){tp30_2<-append(tp30_2,mean(p30_2[,i]))}
for(i in 4:l){twt0_2<-append(twt0_2,mean(wt0_2[,i]))}
for(i in 4:l){twt30_2<-append(twt30_2,mean(wt30_2[,i]))}

#plots
postscript("5ss.eps",height=600,width=600)
#postscript("5ss_with_wt30.eps",height=600,width=600)
plot(tp0,type="l",col=6,main="5'ss",sub=paste("genes=",g1),ylab="reads per base per million mapped reads",ylim=c(min(c(tp0,tp30,twt0,twt30)),max(c(tp0,tp30,twt0,twt30))),axes=F,xlab="Distance to 5'ss")
lines(tp30,col=2)
lines(twt0,col=3)
#lines(twt30,col=5)
axis(1,at=c(1,10,20,30,40,50,60,70,80,90,100),labels=c("-500","-400","-300","-200","-100","0","100","200","300","400","500"))
axis(2,)
legend(1,4.5,legend=c("MT25","MT37","WT"),col=c(6,2,3),lty=c(1,1,1),cex=0.8)
#legend(1,5,legend=c("p0","p30","wt0","wt30"),col=c(2,3,4,5),lty=c(1,1,1,1),cex=0.8)
abline(v=50,col="black")
dev.off()

postscript("5ss_200.eps",height=600,width=600)
#postscript("5ss_200_with_wt30.eps",height=600,width=600)
plot(tp0_2,type="l",col=6,main="5'ss (reads extended to 200b)",sub=paste("genes=",g2),ylab="reads per base per million mapped reads",ylim=c(min(c(tp0_2,tp30_2,twt0_2,twt30_2)),max(c(tp0_2,tp30_2,twt0_2,twt30_2))),axes=F,xlab="Distance to 5'ss")
lines(tp30_2,col=2)
lines(twt0_2,col=3)
#lines(twt30_2,col=5)
axis(1,at=c(1,10,20,30,40,50,60,70,80,90,100),labels=c("-500","-400","-300","-200","-100","0","100","200","300","400","500"))
axis(2,)
legend(1,20,legend=c("MT25","MT37","WT"),col=c(6,2,3),lty=c(1,1,1),cex=0.8)
#legend(1,5,legend=c("p0","p30","wt0","wt30"),col=c(2,3,4,5),lty=c(1,1,1,1),cex=0.8)
abline(v=50,col="black")
dev.off()

postscript("5ss_log2.eps",height=600,width=600)
plot(log2(tp30/tp0),type="l",col=6,main="5'ss",sub=paste("genes=",g1),ylab="log2(rpbpm(A)/rpbpm(B)",ylim=c(-1,1.5),axes=F,xlab="Distance to 5'ss")
lines(log2(tp0/twt0),col=2)
lines(log2(tp30/twt0),col=3)
axis(1,at=c(1,10,20,30,40,50,60,70,80,90,100),labels=c("-500","-400","-300","-200","-100","0","100","200","300","400","500"))
axis(2,)
legend(80,1.5,legend=c("MT37/MT25","MT25/WT","MT37/WT"),col=c(6,2,3),lty=c(1,1,1),cex=0.8)
abline(0,0,col="black")
abline(v=50,col="black")
dev.off()

postscript("5ss_200_log2.eps",height=600,width=600)
plot(log2(tp30_2/tp0_2),type="l",col=6,main="5'ss (reads extended to 200b)",sub=paste("genes=",g2),ylab="log2(rpbpm(A)/rpbpm(B)",ylim=c(-1,1),axes=F,xlab="Distance to 5'ss")
lines(log2(tp0_2/twt0_2),col=2)
lines(log2(tp30_2/twt0_2),col=3)
axis(1,at=c(1,10,20,30,40,50,60,70,80,90,100),labels=c("-500","-400","-300","-200","-100","0","100","200","300","400","500"))
axis(2,)
legend(80,1,legend=c("MT37/MT25","MT25/WT","MT37/WT"),col=c(6,2,3),lty=c(1,1,1),cex=0.8)
abline(0,0,col="black")
abline(v=50,col="black")
dev.off()

postscript("5ss_boxplots.eps",height=1000,width=1000)
par(mfrow=c(3,1))
boxplot(p0[,4:l],ylim=c(0,15),col=6,outline=F,axes=F,ylab="rpbpm",xlab="5'SS +/- 500 bases")
axis(2)
abline(v=50,col="black")
title("MT25 5ss")
boxplot(p30[,4:l],ylim=c(0,15),col=2,outline=F,axes=F,ylab="rpbpm",xlab="5'SS +/- 500 bases")
axis(2)
abline(v=50,col="black")
title("MT37 5ss")
boxplot(wt0[,4:l],ylim=c(0,15),col=3,outline=F,axes=F,ylab="rpbpm",xlab="5'SS +/- 500 bases")
axis(2)
abline(v=50,col="black")
title("WT 5ss")
#boxplot(wt30[,4:l],ylim=c(0,15),col=5,outline=F,axes=F,ylab="rpbpm",xlab="5'SS +/- 500 bases")
#axis(2)
#abline(v=50,col="black")
#title("wt30 5ss")
dev.off()

postscript("5ss_200_boxplots.eps",height=1000,width=1000)
par(mfrow=c(3,1))
boxplot(p0_2[,4:l],ylim=c(0,30),col=6,outline=F,axes=F,ylab="rpbpm",xlab="5'SS +/- 500 bases")
axis(2)
abline(v=50,col="black")
title("MT25 5ss")
boxplot(p30_2[,4:l],ylim=c(0,30),col=2,outline=F,axes=F,ylab="rpbpm",xlab="5'SS +/- 500 bases")
axis(2)
abline(v=50,col="black")
title("MT37 5ss")
boxplot(wt0_2[,4:l],ylim=c(0,30),col=3,outline=F,axes=F,ylab="rpbpm",xlab="5'SS +/- 500 bases")
axis(2)
abline(v=50,col="black")
title("WT 5ss")
#boxplot(wt30_2[,4:l],ylim=c(0,30),col=5,outline=F,axes=F,ylab="rpbpm",xlab="5'ss +/- 500 bases")
#axis(2)
#abline(v=50,col="black")
#title("wt30 5ss")
dev.off()