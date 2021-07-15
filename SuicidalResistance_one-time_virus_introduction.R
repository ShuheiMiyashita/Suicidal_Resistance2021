#Text S4
#An R script for simulation of suicidal population resistance of SHR(+) plants with one-time virus introduction

## CAUTION: Many files will be produced automatically. Before you start simulation, please (make and) select a directory to which the files will be saved. To complete the simulation below, several hours may be required, if you use a standard laptop. A modified R script for testing parameters manually without producing files (but showing Fig. 6-like results) can be found as Text S3.

#### initial settings
s <- 100      # size of lattice
rp <- 0.3     # reproduction rate of plant
rv <- 1.2     # reproduction rate of virus
lv <- 0.9     # dependency on local reproduction of virus
d <- 0.1      # mortality of plant without infection
dni <- 0.2    # mortality of SHR(-) plant upon viral infection 
dpi <- 1      # mortality of SHR(+) plant upon viral infection
vo <- 0.1     # one-time virus introduction


#### simulations in different conditions for Rshrp and lp
for (Rshrp in c(0.1,0.5,0.9)){ #Rshrp indicates initial proportion of SHR(+) plants
pt20rec <- NULL #for recording extinction time
ptcrec <- NULL #for recording cumulative number of SHR(-) and SHR(+) plants 
for (lp in c(1,0.95,0.9)){ #lp indicates dependency on local reproduction of plants
for (seed in 1:10){ #10 trials
###initial settings for each trial
set.seed(seed)
rptv <- NULL # for recording plant and virus abundance at different time points
pt20rec <- c(pt20rec,lp,seed) #registration of current trial
pt20 <- 1 # 1: before extinction; 0: after extinction of SHR(+) plants
ptc <- c(0,0) #for recording cumulative number of SHR(-) and SHR(+) plants

### running-in without virus
tr <- 0 # time after starting running-in
pt <- matrix(rep(1,s*s),nrow=s)  #plant table; 0: open box, 1: SHR(-), and 2: SHR(+)
pt[sample(1:(s*s),round(s*s*Rshrp,0),replace=F)] <- 2  #introducing SHR(+) plants according to lp
while (tr < 201){
pt <- pt*rbinom(s*s,1,1-d) # death without virus
## visualization
if (tr%%50==0){
pname <- paste("plot_o",seed,"_",Rshrp*100,"_",lp*100,"_tr",tr,".png",sep="")
png(pname,width=800,height=800)
plot(0,0,xlim=c(0,s),ylim=c(0,s),type="n",xlab="",ylab="")
for (x in 1:s){
points(rep(x,s),1:s,col=rgb((1-pt[x,])^2,1-pt[x,]*(pt[x,]-1)/2,(1-pt[x,])^2),pch=19,cex=1) # SHR(-): green; SHR(+): magenta
}
dev.off()
ptv <- c(tr,length(which(pt==0)),length(which(pt==1)),length(which(pt==2)),0)
rptv <- rbind(rptv,ptv)
}else{
}
## plant propagation
pto <- (pt-1)*(pt-2)/2                                   # finding open box
plt <- pt[c(s,1:(s-1)),]*rbinom(s*s,1,rp*lp/4)*pto       # local propagation from top box
plb <- pt[c(2:s,1),]*rbinom(s*s,1,rp*lp/4)*pto           # local propagation from bottom box
pll <- pt[,c(s,1:(s-1))]*rbinom(s*s,1,rp*lp/4)*pto       # local propagation from left box
plr <- pt[,c(2:s,1)]*rbinom(s*s,1,rp*lp/4)*pto           # local propagation from right box
pgp <- sample(1:(s*s),rbinom(1,sum((pt-1)*pt/2),rp*(1-lp)),replace=F); pgpv <- rep(0,s*s); pgpv[pgp] <- 2 # global propagation of SHR(+) plant 
pgn <- sample(1:(s*s),rbinom(1,sum((pt-2)^2*pt),rp*(1-lp)),replace=F); pgnv <- rep(0,s*s); pgnv[pgn] <- 1 # global propagation of SHR(-) plant
## deciding which parent leave progeny
pv <- rbind(as.vector(plt),as.vector(plb),as.vector(pll),as.vector(plr),pgpv*as.vector(pto),pgnv*as.vector(pto))
cspv <- colSums(pv)
pvp <-  1-(pv-1)*(pv-2)/2                      # returns 1 if a box has plant inhabitant
cspvp <- colSums(pvp)                          # number of inhabited boxes
prv <- rep(0,s*s)                              # vector for propagation result
prv[which(cspvp==1)] <- cspv[which(cspvp==1)]  # boxes with only one parent candidate
for (i in which(cspvp > 1)){                   # boxes with multiple parent candidates
prv[i] <- sample(pv[,i],1,prob=pvp[,i])        # random decision of parents
}
pt <- pt+prv
tr <- tr+1
}

### main part after starting virus introduction
t <- 0 # time after starting virus introduction
vt <- matrix(rep(0,s*s),nrow=s) # generating virus table and virus introduction; 0: no infection and 1: infection
tvo <- rbinom(s*s,1,vo)*(1-vt)*(1-(pt-1)*(pt-2)/2) #10% of plants are randomly challenged by the virus
vt <- vt+tvo


while (t < 5001){
## death
drt <- ((1-vt)*d+vt*pt*(pt-1)^2*dpi/2+vt*pt*(pt-2)^2*dni)*ceiling(pt/2) # mortality table
pt <- pt*rbinom(s*s,1,1-drt)
## viral propagation
pti <- 1-(pt-1)*(pt-2)/2                                     # boxes inhabited
vt <- vt*pti                                                 # updating vt by removing dead plant
vlt <- vt[c(s,1:(s-1)),]*rbinom(s*s,1,rv*lv/4)               # local spread from top box
vlb <- vt[c(2:s,1),]*rbinom(s*s,1,rv*lv/4)                   # local spread from bottom box
vll <- vt[,c(s,1:(s-1))]*rbinom(s*s,1,rv*lv/4)               # local spread from left box
vlr <- vt[,c(2:s,1)]*rbinom(s*s,1,rv*lv/4)                   # local spread from right box
vg <- sample(1:(s*s),rbinom(1,sum(vt),rv*(1-lv)),replace=F)  # global spread 
vgv <- rep(0,s*s); vgv[vg] <- 1                              # global spread vector
vv <- as.vector(vlt)+as.vector(vlb)+as.vector(vll)+as.vector(vlr)+vgv # sum
vt <- vt+vv                                                  # addition allowing >1
vt <- ceiling(vt/6)*pti     # updating vt by limiting to plant-inhabiting boxes, without allowing >1
## plots and record
if (t%%50==0){
pname <- paste("plot_o",seed,"_",Rshrp*100,"_",lp*100,"_t",t,".png",sep="")
png(pname,width=800,height=800)
plot(0,0,xlim=c(0,s),ylim=c(0,s),type="n",xlab="",ylab="")
for (x in 1:s){
points(rep(x,s),1:s,col=rgb((1-pt[x,])^2,1-pt[x,]*(pt[x,]-1)/2,(1-pt[x,])^2),pch=19,cex=1) # SHR(-): green; SHR(+): magenta
points(rep(x,s),1:s,col=rgb(0,0,0,alpha=vt[x,]),pch=22,cex=1)
}
dev.off()
ptv <- c(200+t,length(which(pt==0)),length(which(pt==1)),length(which(pt==2)),sum(vt))
rptv <- rbind(rptv,ptv)
}else{
}
## plant propagation
pto <- (pt-1)*(pt-2)/2 # finding open box
plt <- pt[c(s,1:(s-1)),]*rbinom(s*s,1,rp*lp/4)*pto       # local propagation from top box
plb <- pt[c(2:s,1),]*rbinom(s*s,1,rp*lp/4)*pto           # local propagation from bottom box
pll <- pt[,c(s,1:(s-1))]*rbinom(s*s,1,rp*lp/4)*pto       # local propagation from left box
plr <- pt[,c(2:s,1)]*rbinom(s*s,1,rp*lp/4)*pto           # local propagation from right box
pgp <- sample(1:(s*s),rbinom(1,sum((pt-1)*pt/2),rp*(1-lp)),replace=F); pgpv <- rep(0,s*s); pgpv[pgp] <- 2 # global propagation of SHR(+) plant 
pgn <- sample(1:(s*s),rbinom(1,sum((pt-2)^2*pt),rp*(1-lp)),replace=F); pgnv <- rep(0,s*s); pgnv[pgn] <- 1 # global propagation of SHR(-) plant
## deciding which parent leave progeny
pv <- rbind(as.vector(plt),as.vector(plb),as.vector(pll),as.vector(plr),pgpv*as.vector(pto),pgnv*as.vector(pto))
cspv <- colSums(pv)
pvp <-  1-(pv-1)*(pv-2)/2                      # returns 1 if a box has plant inhabitant
cspvp <- colSums(pvp)                          # number of inhabited boxes
prv <- rep(0,s*s)                              # vector for propagation result
prv[which(cspvp==1)] <- cspv[which(cspvp==1)]  # boxes with only one parent candidate
for (i in which(cspvp > 1)){                   # boxes with multiple parent candidates
prv[i] <- sample(pv[,i],1,prob=pvp[,i])        # random decision of parent
}
pt <- pt+prv
ptc <- ptc+c(length(which(pt==1)),length(which(pt==2))) # updating cumulative number
## detecting and recording extinction
if(length(which(pt==2))==0&&pt20==1){
pt20rec <- c(pt20rec,t)
pt20 <- 0
}else{
}
t <- t+1
}
ptcrec <- rbind(ptcrec,ptc) #recording cumulative numbers

### outputting the record of plant and virus abundance for each trial
fname <- paste("output_o",seed,"_",Rshrp*100,"_",lp*100,".csv",sep="")
write.csv(rptv,file=fname,row.names=F)

}

}

### outputting the record of extinction time
cname1 <- paste("pt20rec_o_",Rshrp*100,".csv",sep="")
write.csv(pt20rec,file=cname1,row.names=F)

### outputting the record of cumulative number
cname2 <- paste("ptcrec_o_",Rshrp*100,".csv",sep="")
write.csv(ptcrec,file=cname2,row.names=F)
}


#### summary of outputs
###vidualization
for (k in c(10,50,90)){
for (j in c(90,95,100)){
pname <- paste("fig_o",k,"_",j,".png",sep="")
png(pname,width=800,height=800)
plot(0,0,xlim=c(0,105),ylim=c(0,1),type="n",xlab="time",ylab="proportion of SHR(+)",xaxt="n",yaxt="n")
axis(side=1,at=0:21*5,labels=F)
axis(side=2,at=0:10*0.1,labels=F)
axis(side=4,at=0:10*0.1,labels=F)
for (i in 1:10){ 
fname <- paste("output_o",i,"_",k,"_",j,".csv",sep="")
tt <- read.csv(fname)
par(new=T)
plot(0:105,tt[1:106,4]/(tt[1:106,3]+tt[1:106,4]),pch=19,xlim=c(0,105),ylim=c(0,1),col=rgb(1,0,1,alpha=0.3),type="l",lwd=3,ann=F,xlab="",ylab="",xaxt="n",yaxt="n")
par(new=T)
plot(5:105,tt[6:106,5]/(tt[6:106,3]+tt[6:106,4]),pch=19,xlim=c(0,105),ylim=c(0,1),col=rgb(0,0,0,alpha=0.2),type="l",lwd=3,ann=F,xlab="",ylab="",xaxt="n",yaxt="n")
}
abline(v=5,lty=2)
abline(v=0,lty=2)
dev.off()
}
}

### occupancy
tvr <- NULL
for (k in c(10,50,90)){
for (j in c(90,95,100)){
tv <- NULL
for (i in 1:10){ 
fname <- paste("output_o",i,"_",k,"_",j,".csv",sep="")
tt <- read.csv(fname)
tv <- c(tv,tt[106,4]/(tt[106,3]+tt[106,4]))
}
tvr <- rbind(tvr,c(mean(tv),sd(tv)))
}
}
write.csv(tvr,file="final_occupancy_o.csv")

