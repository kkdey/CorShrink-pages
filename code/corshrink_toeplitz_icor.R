library(CorShrink)
library(glasso)
library(corpcor)
library(psych)
library(Matrix)

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


frob_vals <- matrix(0, NUM_SIM, 3)

for(m in 1:NUM_SIM){
  ll <- DM_toeplitz(n=n, P=P)
  data <- rbind(ll$Xtrain, ll$Xtest)
  Sigma <- ll$Sigma
  corSigma <- cov2cor(Sigma)
  icorSigma <- solve(corSigma)

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

  strimmer_sample <- corpcor::invcor.shrink(data)

  ## CorShrink

  cov_sample_ML <-  CorShrinkData(data, sd_boot = FALSE,
                                  type = "cor", cor_method = "pearson",
                                  image = "null",
                                  ash.control = list())
  icor_shrink  <- solve(cov_sample_ML$ash_cor_PD)


  frob_glasso <- mean(as.matrix(cov2cor(a$w)) -as.matrix(cov2cor(corSigma)))^2
  #frob_glasso <- 1 - (tr(as.matrix(a$wi %*% icorSigma)))/(norm(a$wi, type = "F")* norm(icorSigma, type = "F"))

  frob_strimmer <- mean(as.matrix(cov2cor(strimmer_sample)) - as.matrix(cov2cor(corSigma)))^2
  #frob_strimmer <- 1 - (tr(as.matrix(strimmer_sample[1:P,1:P]%*%icorSigma)))/(norm(strimmer_sample[1:P,1:P], type = "F")* norm(icorSigma, type = "F"))

  frob_corshrink <- mean(as.matrix(cov2cor(cov_sample_ML$ash_cor_PD)) - as.matrix(cov2cor(corSigma)))^2
  #frob_corshrink <- 1 - (tr(as.matrix(icor_shrink%*%icorSigma)))/(norm(icor_shrink, type = "F")* norm(icorSigma, type = "F"))

  frob_vals[m, ] <- c(frob_glasso, frob_strimmer, frob_corshrink)
  cat("We are at simulation", m, "\n")
}

save(frob_vals, file = paste0("toeplitz_cmd_n_", n, "_P_", P, "_results_icor.rda"))


