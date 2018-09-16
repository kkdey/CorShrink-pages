

#############  CorShrink on Toeplitz matrices   ######################

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


frob_vals <- matrix(0, NUM_SIM, 4)

for(m in 1:NUM_SIM){
  ll <- DM_toeplitz(n=n, P=P)
  data <- rbind(ll$Xtrain, ll$Xtest)
  Sigma <- ll$Sigma
  corSigma <- cov2cor(Sigma)

  ###  GLASSO

  nr  <- 100
  rho <- seq(1e-04,10,length=nr)
  bic <- rho
  S <- cov(data)
  for(j in 1:nr){
    a       <- glasso::glasso(S,rho[j])
    p_off_d <- sum(a$wi!=0 & col(S)<row(S))
    bic[j]  <- -2*(a$loglik) + p_off_d*log(n)
  }
  best <- which.min(bic)

  a <- glasso::glasso(S,rho[best])

  ## Strimmer Shafer

  strimmer_sample <- corpcor::cor.shrink(data)

  ## CorShrink

  cov_sample_ML <-  CorShrinkData(data, sd_boot = FALSE,
                                  ash.control = list())


  # frob_S <- mean(as.matrix(cov2cor(S)) - as.matrix(cov2cor(corSigma)))^2
  frob_S <- 1 - (tr(as.matrix(cov2cor(S)%*%cov2cor(corSigma))))/(norm(cov2cor(S), type = "F")* norm(corSigma, type = "F"))

  # frob_glasso <- mean(as.matrix(cov2cor(a$w)) -as.matrix(cov2cor(corSigma)))^2
  frob_glasso <- 1 - (tr(as.matrix(cov2cor(a$w)%*%cov2cor(corSigma))))/(norm(cov2cor(a$w), type = "F")* norm(corSigma, type = "F"))

  # frob_strimmer <- mean(as.matrix(cov2cor(strimmer_sample)) - as.matrix(cov2cor(corSigma)))^2
  frob_strimmer <- 1 - (tr(as.matrix(cov2cor(strimmer_sample[1:P,1:P])%*%cov2cor(corSigma))))/(norm(cov2cor(strimmer_sample[1:P,1:P]), type = "F")* norm(corSigma, type = "F"))

  # frob_corshrink <- mean(as.matrix(cov2cor(cov_sample_ML$ash_cor_PD)) - as.matrix(cov2cor(corSigma)))^2
  frob_corshrink <- 1 - (tr(as.matrix(cov2cor(cov_sample_ML$ash_cor_PD)%*%cov2cor(corSigma))))/(norm(cov2cor(cov_sample_ML$ash_cor_PD), type = "F")* norm(corSigma, type = "F"))


  frob_vals[m, ] <- c(frob_S, frob_glasso, frob_strimmer, frob_corshrink)
  cat("We are at simulation", m, "\n")
}

save(frob_vals, file = paste0("toeplitz_cmd_n_", n, "_P_", P, "_results.rda"))


