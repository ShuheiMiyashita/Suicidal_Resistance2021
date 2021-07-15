#Text S3
#A modified R script for simulation of suicidal population resistance of SHR(+) plants, without producing files

## This script is for testing different parameters manually. For making multiple trials with different lp or Rshrp values, please refer Text S3. 

#### initial settings
s <- 100      # size of lattice
Rshrp <- 0.5  # initial proportion of SHR(+) plants
rp <- 0.3     # reproduction rate of plant
rv <- 1.2     # reproduction rate of virus
lp <- 1.0     # dependency on local reproduction of plants
lv <- 0.9     # dependency on local reproduction of virus
d <- 0.1      # mortality of plant without infection
dni <- 0.2    # mortality of SHR(-) plant upon viral infection 
dpi <- 1      # mortality of SHR(+) plant upon viral infection
vi <- 0.001   # rate of virus influx


#### main body of simulation
###initial settings for each trial
rptv <- NULL # for recording plant and virus abundance at different time points

### running-in without virus
tr <- 0 # time after starting running-in
pt <- matrix(rep(1,s*s),nrow=s)  #plant table; 0: open box, 1: SHR(-), and 2: SHR(+)
pt[sample(1:(s*s),round(s*s*Rshrp,0),replace=F)] <- 2  #introducing SHR(+) plants according to lp
while (tr < 201){
pt <- pt*rbinom(s*s,1,1-d) # death without virus
## visualization
if (tr%%50==0){
plot(0,0,xlim=c(0,s),ylim=c(0,s),type="n",xlab="",ylab="")
for (x in 1:s){
points(rep(x,s),1:s,col=rgb((1-pt[x,])^2,1-pt[x,]*(pt[x,]-1)/2,(1-pt[x,])^2),pch=19,cex=0.5) # SHR(-): green; SHR(+): magenta
}
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

while (t < 5001){
##constant virus influx
cvi <- rbinom(s*s,1,vi)*(1-vt)*(1-(pt-1)*(pt-2)/2) #0.1% of plants are randomly challenged by the virus
vt <- vt+cvi
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
plot(0,0,xlim=c(0,s),ylim=c(0,s),type="n",xlab="",ylab="")
for (x in 1:s){
points(rep(x,s),1:s,col=rgb((1-pt[x,])^2,1-pt[x,]*(pt[x,]-1)/2,(1-pt[x,])^2),pch=19,cex=0.5) # SHR(-): green; SHR(+): magenta
points(rep(x,s),1:s,col=rgb(0,0,0,alpha=vt[x,]),pch=22,cex=0.5)
}
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
t <- t+1
}
