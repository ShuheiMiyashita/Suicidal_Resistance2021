#Text S1
#An R script for MOI estimation

### Data input
m0 <- 63; s0 <- 7; sites012 <- 28; m012 <- 457; s012 <- 172    #for N. benthamiana(R-) - wildtype RNA3
#NOTE: select one line for other combinations by removing first "#" in each line:
#m0 <- 27; s0 <- 6; sites012 <- 10; m012 <- 197; s012 <- 78    #for N. benthamiana(R-) - RNA3 CP-N31T
#m0 <- 29; s0 <- 9; sites012 <- 11; m012 <- 214; s012 <- 84    #for N. benthamiana(R-) - RNA3 CP-T45M 
#m0 <- 40; s0 <- 27; sites012 <- 19; m012 <- 232; s012 <- 156  #for N. benthamiana(R+) - wildtype RNA3 
#m0 <- 45; s0 <- 16; sites012 <- 18; m012 <- 322; s012 <- 125  #for N. benthamiana(R+) - RNA3 CP-N31T
#m0 <- 36; s0 <- 17; sites012 <- 15; m012 <- 258; s012 <- 130  #for N. benthamiana(R+) - RNA3 CP-T45M 

# m0: no. of co-infected sites
# s0: no. of singly infected sites 
# sites12: no. of co-infected sites with 10-30 infected cells at 14-16 hours after inoculation, used for lambda12 estimation
# m012: no. of co-infected cells
# s012: no. of singly infected cells
n012 <- m012+s012

### Main body of MOI estimation
n012 <- m012+s012
cell0 <- sites012               #number of cell-0 cells in sites012 
cell1 <- sites012*8             #number of cell-1 cells in sites012
cell2 <- n012-sites012*9   #number of cell-2 cells in sites012

K <- 30 #maximum number of founders included in the calculation 
kv <- NULL
lv <- NULL
for (i in 0:K){
kv <- c(kv,rep(i,i+1))
lv <- c(lv,0:i)
}
klv <- kv-lv
lklv <- lv*klv
ln <- (K+2)*(K+1)/2

## Function for log likelihood for any la,bda0 and lambda12
MOImLL <- function(lambda){
lambda0 <- lambda[1]; lambda1 <- lambda[2]; lambda2 <- lambda[2]; r0 <- 0.5
table0 <- matrix(rep(0,ln*5),ncol=5)
table0[,1] <- kv
table0[,2] <- lv
table0[,3] <- klv
table0[,4] <- lklv
table0[,5] <- dpois(table0[,1],lambda0)*dbinom(table0[,2],table0[,1],r0)
pni0 <- table0[1,5]
py0 <- sum(table0[which(table0[,2]>0&table0[,3]==0),5])/(1-pni0)
pc0 <- sum(table0[which(table0[,2]==0&table0[,3]>0),5])/(1-pni0)
pm0 <- 1-py0-pc0
ps0 <- py0+pc0
table0m <- table0[which(table0[,4]>0),]
table0m[,5] <- table0m[,5]/sum(table0m[,5])
table1 <- matrix(rep(0,ln*5),ncol=5)
table1[,1] <- kv
table1[,2] <- lv
table1[,3] <- klv
table1[,4] <- lklv
for (l in 1:nrow(table0m)){
table1[,5] <- table1[,5]+table0m[l,5]*dpois(table1[,1],lambda1)*dbinom(table1[,2],table1[,1],table0m[l,2]/table0m[l,1])
}
pni1 <- table1[1,5]
py1 <- sum(table1[which(table1[,2]>0&table1[,3]==0),5])/(1-pni1)
pc1 <- sum(table1[which(table1[,2]==0&table1[,3]>0),5])/(1-pni1)
pm1 <- 1-py1-pc1
ps1 <- py1 + pc1
table1i <- table1[2:ln,]
table1i[,5] <- table1i[,5]/sum(table1i[,5])
table2 <- matrix(rep(0,ln*5),ncol=5)
table2[,1] <- kv
table2[,2] <- lv
table2[,3] <- klv
table2[,4] <- lklv
for (l in 1:(ln-1)){
table2[,5] <- table2[,5]+table1i[l,5]*dpois(table2[,1],lambda2)*dbinom(table2[,2],table2[,1],table1i[l,2]/table1i[l,1])
}
pni2 <- table2[1,5]
py2 <- sum(table2[which(table2[,2]>0&table2[,3]==0),5])/(1-pni2)
pc2 <- sum(table2[which(table2[,2]==0&table2[,3]>0),5])/(1-pni2)
pm2 <- 1-py2-pc2
ps2 <- py2 + pc2

pm012 <- (cell0+cell1*pm1+cell2*pm2)/n012
ps012 <- (cell1*ps1+cell2*ps2)/n012
logL0 <- dmultinom(c(s0,m0),prob=c(ps0,pm0),log=TRUE)
logL012 <- dmultinom(c(s012,m012),prob=c(ps012,pm012),log=TRUE)
LL <- logL0 + logL012
-LL
}

## Maximization of log likelihood to find most likely lambda
init <- c(5,5)
MOImLL.opt <- optim(init,MOImLL, NULL, method="L-BFGS-B", hessian = TRUE, lower=c(0.1,0.1), upper=c(10,10))
MOImLL.opt$par  # lambda0 and lambda12 estimates
v <- solve(MOImLL.opt$hessian)
se <- sqrt(diag(v))
se              # lambda0 and lambda12 standard errors
