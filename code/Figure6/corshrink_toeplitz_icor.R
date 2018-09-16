library(CorShrink)
library(glasso)
library(corpcor)
library(psych)
library(Matrix)

n <- 1000
P <- 100
NUM_SIM <- 50


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


frob_vals <- matrix(0, NUM_SIM, 7)

for(m in 1:NUM_SIM){
  ll <- DM_toeplitz(n=n, P=P)
  data <- rbind(ll$Xtrain, ll$Xtest)
  Sigma <- ll$Sigma
  corSigma <- cov2cor(Sigma)
  icorSigma <- solve(corSigma)

  S <- cov(data)

  ###  GLASSO

  nr  <- 100
  rho <- c(0.01, 0.1, 0.5, 1)
  a_list <- list()
  for(i in 1:length(rho)){
    a_list[[i]] <- glasso::glasso(S,rho[i])
  }


  ## Strimmer Shafer

  strimmer_sample <- corpcor::invcor.shrink(data)

  ## CorShrink

  cov_sample_ML <-  CorShrinkData(data, sd_boot = FALSE,
                                  ash.control = list())
  icor_shrink  <- solve(cov_sample_ML$cor)

  pdsoft_sample <- PDSCE::pdsoft.cv(data, tolin = 1e-04, tolout = 1e-04)
  ipdsoft <- cov2cor(solve(pdsoft_sample$sigma))


  frob_strimmer <-   mean(sqrt((as.matrix(cov2cor(strimmer_sample[1:P,1:P])) - as.matrix(cov2cor(icorSigma)))^2))

  frob_corshrink <- mean(sqrt((as.matrix(cov2cor(icor_shrink)) - as.matrix(cov2cor(icorSigma)))^2))

  frob_pdsoft <- mean(sqrt((as.matrix(cov2cor(ipdsoft)) - as.matrix(cov2cor(icorSigma)))^2))

  frob_glasso_1 <- mean(sqrt((as.matrix(cov2cor(a_list[[1]]$wi)) - as.matrix(cov2cor(icorSigma)))^2))
  frob_glasso_2 <- mean(sqrt((as.matrix(cov2cor(a_list[[2]]$wi)) - as.matrix(cov2cor(icorSigma)))^2))
  frob_glasso_3 <- mean(sqrt((as.matrix(cov2cor(a_list[[3]]$wi)) - as.matrix(cov2cor(icorSigma)))^2))
  frob_glasso_4 <- mean(sqrt((as.matrix(cov2cor(a_list[[4]]$wi)) - as.matrix(cov2cor(icorSigma)))^2))
 
  frob_vals[m, ] <- c(frob_strimmer, frob_corshrink,  frob_pdsoft,
                      frob_glasso_1,
                      frob_glasso_2,  frob_glasso_3,  frob_glasso_4)
  cat("We are at simulation", m, "\n")
}

save(frob_vals, file = paste0("toeplitz_cmd_n_", n, "_P_", P, "_results_icor.rda"))


