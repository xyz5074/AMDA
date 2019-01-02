############## R function AMDA
## Input: 
## X: a n*p microbiome compositional matrix, each row is a sample (sum to 100%)
## g: a group label of 1 or 2
## B: permutation times 
## Output:
## pval: AMDA test p-value 		 
AMDA <-function(X, g2, B=1000){

getLogCLR = function(comm){
  log_geo_mean <- function(data) {
    log_data <- log(data)
    log_gm <- mean(log_data[is.finite(log_data)])  # log 0 issue..but in simulation we have added 0.01 to each value
    return(log_gm)
  }
  log_geoMean = apply(comm, 1, log_geo_mean)
  logclr= log(comm) - log_geoMean
  return(logclr)  
}

kernel.gauss<-function(zz){
## Input: zz: a n*p data matrix, each row is a sample, each column is a variable
## Output:Kmat: a n*n kernel matrix
	if (is.null(nrow(zz))) {zz=as.matrix(zz,ncol=1)} ##zz is a 1-row vector
	n<-nrow(zz)
	Dmat<-matrix(NA, nrow=n, ncol=n)
	for (i in 1:n){
		for(j in 1:n){
			Dmat[i,j]<- sum((zz[i,]-zz[j,])^2)
		}
	}
      scl=median(Dmat)
	Kmat<-matrix(NA, nrow=n, ncol=n)
	for (i in 1:n){
		for(j in 1:n){
			Kmat[i,j]<-exp(-sum((zz[i,]-zz[j,])^2)/scl)
		}
	}	
	return(Kmat)		
}

############## MMD   
## Input: 
## K: a n*n kernel matrix
## g: a group label of 1 or 2
## B: permutation times 
## Output:
## stat: MMD test statistic
## pval: MMD test p-value (only if Test=TRUE)		 
MMD<-function(K, g, Test=FALSE, B=1000){
n=nrow(K)
K11=K[g==1,g==1]
K12=K[g==1,g==2]
K21=K[g==2,g==1]
K22=K[g==2,g==2]
MMD=mean(K11)+mean(K22)-mean(K12)-mean(K21)
if (!Test) { return(MMD) }
if (Test){
MMDperm=rep(NA,B)
for(j in 1:B){
pindex=sample(1:n,size=n,replace=F)
gp=g[pindex]
K11=K[gp==1,gp==1]
K12=K[gp==1,gp==2]
K21=K[gp==2,gp==1]
K22=K[gp==2,gp==2]
MMDperm[j]=mean(K11)+mean(K22)-mean(K12)-mean(K21)
}
pv=sum(MMDperm>MMD)/B
return(list(stat=MMD,pval=pv))
}}

n=nrow(X)
p=ncol(X)
Z = getLogCLR(X)
Z_p=Z[c(sample(nrow(Z))), ]
pvec=pvecp=rep(NA,p)
for(j in 1:p){ pvec[j]=ks.test(Z[g2==1,j],Z[g2==2,j])$p.value}
for(j in 1:p){ pvecp[j]=ks.test(Z_p[g2==1,j],Z_p[g2==2,j])$p.value}
subt=NULL
for(i in 1:p){if(pvec[i] < pvecp[i]) subt=c(subt,i)}
if(length(subt)==0){ subt = order(pvec)[1]} ## 0 selected is possible for p=20
kgsub=kernel.gauss(Z[,subt])
Tobs=MMD(kgsub,g2,Test=FALSE)
T=rep(NA, B)
for(t in 1:B){
pindex=sample(1:n,size=n,replace=F)
g2p=g2[pindex]
pvec=pvecp=rep(NA,p)
for(j in 1:p){ pvec[j]=ks.test(Z[g2p==1,j],Z[g2p==2,j])$p.value}
for(j in 1:p){ pvecp[j]=ks.test(Z_p[g2p==1,j],Z_p[g2p==2,j])$p.value}
subt=NULL
for(i in 1:p){if(pvec[i] < pvecp[i] ) subt=c(subt,i)}
if(length(subt)==0){ subt = order(pvec)[1]}
kgsub=kernel.gauss(Z[,subt])
T[t]=MMD(kgsub,g2p,Test=FALSE)
}
pval=sum(T>Tobs)/B
return(pval)
}

## Testing example:
library(MASS)
n1=n2=25
n=n1+n2
g2=c(rep(1,n1),rep(2,n2))
p=20 
mu1=mu2=runif(p,min=0,max=10)
D=diag(sqrt(runif(p,min=1,max=3)),p)
A = matrix(0,p,p)
for(i in 1:p){
for(j in 1:p){
if (i==j) {A[i,j]=1}
if(abs(i-j)==1) {A[i,j]=-0.5}
}}
sigma=D%*%A%*%D
sigma2=sigma1=sigma
W= mvrnorm(n1, mu1, sigma1, tol = 1e-6)
W= exp(W)
rowW = rowSums(W)
colW = colSums(W)
X1 = W
for(i in 1:nrow(W)){
for(j in 1:ncol(W)){
X1[i,j] = W[i,j]/rowW[i]
}}
W= mvrnorm(n2, mu2, sigma2, tol = 1e-6)
W= exp(W)
rowW = rowSums(W)
colW = colSums(W)
X2 = W
for(i in 1:nrow(W)){
for(j in 1:ncol(W)){
X2[i,j] = W[i,j]/rowW[i]
}}
X=rbind(X1,X2)
AMDA(X, g2, B=100)

