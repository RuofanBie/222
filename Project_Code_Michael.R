data=matrix(1:96, nrow=24)
colnames(data)=paste(c("a","Wald","Fisher","Score"))
data=data.frame(data)
data$a=1:24
data$a=data$a*2

#Wald
for (j in 1:24){
  table=matrix(c(data[j,1],50-data[j,1],data[j,1]/2,50-data[j,1]/2),ncol=2,byrow=TRUE)
  colnames(table)=c("Yes","No")
  rownames(table)=c("A","B")
  table=as.table(table)
  wald=oddsratio.wald(table, y=NULL,
                    conf.level=0.95,
                    rev=c("neither"),
                    correction=FALSE,
                    verbose=FALSE)
  wald=wald$measure
  data[j,2]=wald[2,3]-wald[2,2]
}

#Fisher
for (j in 1:24){
  table=matrix(c(data[j,1],50-data[j,1],data[j,1]/2,50-data[j,1]/2),ncol=2,byrow=TRUE)
  colnames(table)=c("Yes","No")
  rownames(table)=c("A","B")
  table=as.table(table)
  fisher=oddsratio.fisher(table, y=NULL,
                      conf.level=0.95,
                      rev=c("neither"),
                      correction=FALSE,
                      verbose=FALSE)
  fisher=fisher$measure
  data[j,3]=fisher[2,3]-fisher[2,2]
}

#Score
for (j in 1:24){
  score=orscoreci(data[j,1],50,data[j,1]/2,50,0.95)
  data[j,4]=score[[1]][2]-score[[1]][1]
}

#reshape to long
data.1=melt(data, id="a")
colnames(data.1)=paste(c("a","Type","Width"))

#graph the three types of CIs
ggplot(data.1, aes(x=a, y=Width, group=Type)) + 
  labs(y="Confidence Interval Width",x="2x2 table") +
  geom_line(aes(color=Type), size=1) +
  scale_x_continuous(breaks=c(10,20,30,40), labels=c("10:40/5:45","20:30/10:40","30:20/15:35","40:10/20:30"))

#a closer look
ggplot(data.1, aes(x=a, y=Width, group=Type)) + 
  coord_cartesian(ylim=c(0,20)) +
  labs(y="Confidence Interval Width",x="2x2 table") +
  geom_line(aes(color=Type), size=1) +
  scale_x_continuous(breaks=c(10,20,30,40), labels=c("10:40/5:45","20:30/10:40","30:20/15:35","40:10/20:30"))

#Michael - Confidence Interval Width Simulations
#Based on hypergeometric, varying the case:control ratio from 2:98 to 98:2, fixing the row totals - 100 iterations each for Wald, Fisher and Score

if (!require("epitools")) install.packages("epitools")
if (!require("PropCIs")) install.packages("PropCIs")
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("reshape")) install.packages("reshape")
library("epitools")
library("PropCIs")
library("ggplot2")
library("reshape")

#create empty matrix and generate cell 'a' for 2x2 based on hypergeometric distribution fixing row totals
data=matrix(1:4900, nrow=100)
colnames(data)=paste("sim",1:49,sep="")
data=data.frame(data)
for (i in 1:49){
  set.seed(i)
  data[1:100,i]=rhyper(100,i*2,100-i*2,50)
}

#create empty df for summary measures to be filled in
summary=matrix(1:196, nrow=49)
colnames(summary)=paste(c("sim","Wald","Fisher","Score"))
summary=data.frame(summary)
summary$sim=1:49

#Loop 49 different scenarios (2:98 to 98:2) for Wald, Fisher, Score - 100 iterations each
for (k in 1:49){

### Wald
sim.wald=matrix(1:200, nrow=100)
colnames(sim.wald)=paste(c("lci","uci"))
for (j in 1:100){
  table=matrix(c(data[j,k],50-data[j,k],k*2-data[j,k],100-k*2-50+data[j,k]),ncol=2,byrow=TRUE)
  colnames(table)=c("Yes","No")
  rownames(table)=c("A","B")
  table=as.table(table)
  wald=oddsratio.wald(table, y=NULL,
                      conf.level=0.95,
                      rev=c("neither"),
                      correction=FALSE,
                      verbose=FALSE)
  wald=wald$measure
  sim.wald[j,1]=wald[2,2]
  sim.wald[j,2]=wald[2,3]
}
sim.wald=data.frame(sim.wald)
sim.wald$width=sim.wald$uci-sim.wald$lci
summary[k,2]=median(sim.wald$width,na.rm=T)

### Fisher
sim.fisher=matrix(1:200, nrow=100)
colnames(sim.fisher)=paste(c("lci","uci"))
for (j in 1:100){
  table=matrix(c(data[j,k],50-data[j,k],k*2-data[j,k],100-k*2-50+data[j,k]),ncol=2,byrow=TRUE)
  colnames(table)=c("Yes","No")
  rownames(table)=c("A","B")
  table=as.table(table)
  fisher=oddsratio.fisher(table, y=NULL,
                      conf.level=0.95,
                      rev=c("neither"),
                      correction=FALSE,
                      verbose=FALSE)
  fisher=fisher$measure
  sim.fisher[j,1]=fisher[2,2]
  sim.fisher[j,2]=fisher[2,3]
}
sim.fisher=data.frame(sim.fisher)
sim.fisher$width=sim.fisher$uci-sim.fisher$lci
summary[k,3]=median(sim.fisher$width,na.rm=T)

### Score
sim.score=matrix(1:200, nrow=100)
colnames(sim.score)=paste(c("lci","uci"))
for (j in 1:100){
  score=orscoreci(data[j,k],50,k*2-data[j,k],50,0.95)
  sim.score[j,1]=score[[1]][1]
  sim.score[j,2]=score[[1]][2]
}
sim.score=data.frame(sim.score)
sim.score$width=sim.score$uci-sim.score$lci
summary[k,4]=median(sim.score$width,na.rm=T)

}

#reshape to long format
summary.1=melt(summary, id="sim")
colnames(summary.1)=paste(c("Sim","Type","Width"))

#graph the three types of CIs
ggplot(summary.1, aes(x=Sim, y=Width, group=Type)) + 
  labs(y="Confidence Interval Width",x="Cases:Controls") +
  geom_line(aes(color=Type), size=1) +
  scale_x_continuous(breaks=c(5,15,25,35,45), labels=c("10:90","30:70","50:50","70:30","90:10"))

#graph but restrict to max y=20
ggplot(summary.1, aes(x=Sim, y=Width, group=Type)) + 
  coord_cartesian(ylim=c(0,25)) +
  labs(y="Confidence Interval Width",x="Cases:Controls") +
  geom_line(aes(color=Type), size=1) +
  scale_x_continuous(breaks=c(5,15,25,35,45), labels=c("10:90","30:70","50:50","70:30","90:10"))
