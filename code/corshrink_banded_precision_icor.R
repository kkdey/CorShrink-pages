
library(CorShrink)
library(corpcor)
library(glasso)
library(Matrix)
library(psych)

n <- 50
P <- 100
NUM_SIM <- 50

diags <- list()
diags[[1]] <- rep(1, 100)
diags[[2]] <- rep(-0.5, 100)
Kinv <- bandSparse(100, k = -(0:1), diag = diags, symm = TRUE)
K <- solve(Kinv)
corSigma <- cov2cor(K)
icorSigma <- solve(corSigma)

frob_vals <- matrix(0, NUM_SIM, 4)

for(m in 1:NUM_SIM){

  data <- MASS::mvrnorm(n,rep(0,P),corSigma)
  S <- cov(data)

  pdsoft_sample <- PDSCE::pdsoft.cv(data, tolin = 1e-04, tolout = 1e-04)
  ipdsoft <- cov2cor(solve(pdsoft_sample$sigma))
  mean(sqrt((as.matrix(cov2cor(ipdsoft)) - as.matrix(cov2cor(icorSigma)))^2))

  ###  GLASSO

  nr  <- 100
  rho <- seq(1e-10,10,length=nr)
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

  frob_glasso <- mean((as.matrix(cov2cor(a$w)) -as.matrix(cov2cor(corSigma)))^2)
  #frob_glasso <- 1 - (tr(as.matrix(a$wi %*% icorSigma)))/(norm(a$wi, type = "F")* norm(icorSigma, type = "F"))

  frob_strimmer <- mean((as.matrix(cov2cor(strimmer_sample)) - as.matrix(cov2cor(corSigma)))^2)
  #frob_strimmer <- 1 - (tr(as.matrix(strimmer_sample[1:P,1:P]%*%icorSigma)))/(norm(strimmer_sample[1:P,1:P], type = "F")* norm(icorSigma, type = "F"))

  frob_corshrink <- mean((as.matrix(cov2cor(cov_sample_ML$ash_cor_PD)) - as.matrix(cov2cor(corSigma)))^2)
  #frob_corshrink <- 1 - (tr(as.matrix(icor_shrink%*%icorSigma)))/(norm(icor_shrink, type = "F")* norm(icorSigma, type = "F"))

  frob_vals[m, ] <- c(frob_glasso, frob_strimmer, frob_corshrink)
  cat("We are at simulation", m, "\n")


}

save(frob_vals, file = paste0("banded_precision_n_", n, "_P_", P, "_results_icor.rda"))
