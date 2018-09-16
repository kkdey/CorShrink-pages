library(CorShrink)
library(glasso)
library(corpcor)
library(Matrix)
library(psych)

n <- 10
P <- 100
NUM_SIM <- 5

block <- 10
mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
Sigma <-   bdiag(mat, mat, mat, mat, mat, mat, mat, mat, mat, mat)
corSigma <- cov2cor(Sigma)
icorSigma <- cov2cor(solve(corSigma))

frob_vals <- matrix(0, NUM_SIM, 3)

for(m in 1:NUM_SIM){

  data <- MASS::mvrnorm(n,rep(0,P),Sigma)
  S <- cov(data, method = "pearson")

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
  iglasso <- cov2cor(a$wi)

  ## Strimmer Shafer

  strimmer_sample <- cov2cor(corpcor::invcor.shrink(data))

  ## CorShrink

  cov_sample_ML <-  CorShrinkData(data, sd_boot = FALSE,
                                  type = "cor", cor_method = "pearson",
                                  image = "null",
                                  ash.control = list())
  icor_shrink  <- cov2cor(solve(cov_sample_ML$cor))




  # frob_S <- mean(as.matrix(cov2cor(S)) - as.matrix(cov2cor(corSigma)))^2
  #  frob_S <- 1 - (tr(as.matrix(cor2pcor(cov2cor(S))%*%pcorSigma)))/(norm(cor2pcor(cov2cor(S)), type = "F")* norm(pcorSigma, type = "F"))

  #frob_glasso <- mean(as.matrix(cov2cor(a$w)) -as.matrix(cov2cor(corSigma)))^2
  frob_glasso <- 1 - (tr(as.matrix(iglasso %*% icorSigma)))/(norm(iglasso, type = "F")* norm(icorSigma, type = "F"))

  #frob_strimmer <- mean(as.matrix(cov2cor(strimmer_sample)) - as.matrix(cov2cor(corSigma)))^2
  frob_strimmer <- 1 - (tr(as.matrix(strimmer_sample[1:P,1:P]%*%icorSigma)))/(norm(strimmer_sample[1:P,1:P], type = "F")* norm(icorSigma, type = "F"))

  #frob_corshrink <- mean(as.matrix(cov2cor(cov_sample_ML$cor)) - as.matrix(cov2cor(corSigma)))^2
  frob_corshrink <- 1 - (tr(as.matrix(icor_shrink%*%icorSigma)))/(norm(icor_shrink, type = "F")* norm(icorSigma, type = "F"))

  frob_vals[m, ] <- c(frob_glasso, frob_strimmer, frob_corshrink)
  cat("We are at simulation", m, "\n")

}

save(frob_vals, file = paste0("hub_cmd_boot_n_", n, "_P_", P, "_results_icor.rda"))

