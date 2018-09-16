

#############  CorShrink on hub matrix  ###################

library(CorShrink)
library(glasso)
library(corpcor)
library(Matrix)
library(psych)

n <- 10
P <- 100
NUM_SIM <- 50

block <- 10
mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
Sigma <-   bdiag(mat, mat, mat, mat, mat, mat, mat, mat, mat, mat)
corSigma <- cov2cor(K)


frob_vals <- matrix(0, NUM_SIM, 4)

for(m in 1:NUM_SIM){

  data <- MASS::mvrnorm(n,rep(0,P),Sigma)
  S <- cov(data, method = "spearman")

  ###  GLASSO

  nr  <- 100
  rho <- seq(1e-04,10,length=nr)
  bic <- rho
  S <- cov(data, method = "spearman")
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

  mean(sqrt((as.matrix(cov2cor(S)) - as.matrix(cov2cor(corSigma)))^2))

  mean(as.matrix(cov2cor(strimmer_sample)) - as.matrix(cov2cor(corSigma)))^2

  mean(as.matrix(cov2cor(a$wi)) -as.matrix(cov2cor(corSigma)))^2

 # frob_S <- mean(sqrt((as.matrix(cov2cor(S)) - as.matrix(cov2cor(corSigma)))^2))
  frob_S <- 1 - (tr(as.matrix(cov2cor(S)%*%cov2cor(corSigma))))/(norm(cov2cor(S), type = "F")* norm(corSigma, type = "F"))

 # frob_glasso <- mean(as.matrix(cov2cor(a$w)) -as.matrix(cov2cor(corSigma)))^2
  frob_glasso <- 1 - (tr(as.matrix(cov2cor(a$w)%*%cov2cor(corSigma))))/(norm(cov2cor(a$w), type = "F")* norm(corSigma, type = "F"))

 # frob_strimmer <- mean(as.matrix(cov2cor(strimmer_sample)) - as.matrix(cov2cor(corSigma)))^2
  frob_strimmer <- 1 - (tr(as.matrix(cov2cor(strimmer_sample[1:P,1:P])%*%cov2cor(corSigma))))/(norm(cov2cor(strimmer_sample[1:P,1:P]), type = "F")* norm(corSigma, type = "F"))

 # frob_corshrink <- mean(as.matrix(cov2cor(cov_sample_ML$ash_cor_PD)) - as.matrix(cov2cor(corSigma)))^2
  frob_corshrink <- 1 - (tr(as.matrix(cov2cor(cov_sample_ML$cor)%*%cov2cor(corSigma))))/(norm(cov2cor(cov_sample_ML$cor), type = "F")* norm(corSigma, type = "F"))


  frob_vals[m, ] <- c(frob_S, frob_glasso, frob_strimmer, frob_corshrink)
  cat("We are at simulation", m, "\n")

}

save(frob_vals, file = paste0("hub_cmd_boot_n_", n, "_P_", P, "_results.rda"))

