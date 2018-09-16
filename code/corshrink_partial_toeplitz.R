
######  corshrink partial correlation hub  ###################


library(CorShrink)
library(huge)
library(corpcor)

n <- 10
P <- 100
NUM_SIM <- 100


DM_toeplitz = function(n,P){
  library("MASS")
  index1=sort(sample(seq(1:n),(n/2)))
  index2=seq(1:n)[-index1]
  Sigmatp=function(P){
    a=array(0,dim=c(P,P))
    for(i in 1:P){
      for(j in 1:P){
        a[i,j]=max(1-0.1*(abs(i-j)),0)
      }
    }
    return(a)
  }
  Sigma = Sigmatp(P)
  data = mvrnorm(n,rep(0,P),Sigma)
  Xtest = data[index2,]
  Xtrain = data[index1,]
  Omega = solve(Sigma)
  return(list(Xtrain = Xtrain, Xtest = Xtest, Sigma = Sigma))
}

ll <- DM_toeplitz(n=n, P=P)
data <- rbind(ll$Xtrain, ll$Xtest)
Sigma <- ll$Sigma
corSigma <- cov2cor(Sigma)

col=c(rep(rgb(0,1,0), 10000),rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)), rep(rgb(1,0,0),10000))

image(as.matrix(corSigma),
      col=col, main=paste0("pop corr matrix"), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))

pcorSigma <- cor2pcor(Sigma)

col=c(rep(rgb(0,1,0), 10000),rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)), rep(rgb(1,0,0),10000))

image(as.matrix(pcorSigma),
      col=col, main=paste0("pop corr matrix"), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))
