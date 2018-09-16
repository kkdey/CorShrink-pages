

library(CorShrink)
library(corpcor)

n <- 30
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
  data = mvrnorm(n,rep(100,P),Sigma)
  Xtest = data[index2,]
  Xtrain = data[index1,]
  Omega = solve(Sigma)
  return(list(Xtrain = Xtrain, Xtest = Xtest, Sigma = Sigma))
}


ll <- DM_toeplitz(n=n, P=P)
data <- rbind(ll$Xtrain, ll$Xtest)
Sigma <- ll$Sigma

sample_cor <- cor(data)

Y <- apply(data, c(1,2) , function(x) rpois(1,exp(x)))
Y <- apply(data, c(1,2) , function(x) rpois(1, x))


pearson_cor <- cor(log(Y+1), method = "pearson")
spearman_cor <- cor(log(Y+1), method = "spearman")

plot(as.vector(Sigma[lower.tri(Sigma)]), as.vector(spearman_cor[lower.tri(spearman_cor)]))
plot(as.vector(Sigma[lower.tri(Sigma)]), as.vector(pearson_cor[lower.tri(pearson_cor)]))

library(ggcorrplot)
ggcorrplot(Sigma)
ggcorrplot(pearson_cor)
ggcorrplot(spearman_cor)


corshrink <- CorShrink::CorShrinkMatrix(pearson_cor, nsamp = 30)
ggcorrplot(corshrink$ash_cor_PD)

corshrink2 <- CorShrink::CorShrinkData(data)
ggcorrplot(corshrink2$ash_cor_PD)

data2 <- huge.npn(log(Y+1))
corshrink3 <- CorShrink::CorShrinkData(data2)
ggcorrplot(corshrink3$ash_cor_PD)

sigma_est <- apply(data, 2, function(x) {
                 y <- var(x) - mean(x)
                 if(y > 0){
                   return(y)
                 }else{
                   return(0)
                 }})

cor_4 <- matrix(0, nrow(Sigma), ncol(Sigma))

 cov_y <- cov(Y)
 for(i in 1:nrow(cov_y)){
   cor_4[i,i] <-  cov_y[i,i]/(sigma_est[i]*sigma_est[i])
   if(i > 1){
     for(j in 1:(i-1)){
       cor_4[i,j] <- cov_y[i,j]/(sigma_est[i]*sigma_est[j])
       cor_4[j,i] <- cor_4[i,j]
     }
   }
 }


ggcorrplot(cov2cor(cor_4))

corshrink4 <- CorShrink::CorShrinkMatrix(cov2cor(cor_4), nsamp = 30)
ggcorrplot(corshrink4$ash_cor_PD)

ggcorrplot(corshrink2$ash_cor_PD)

ggcorrplot(corshrink$ash_cor_PD)

