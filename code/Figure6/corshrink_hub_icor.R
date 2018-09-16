library(CorShrink)
library(glasso)
library(corpcor)
library(Matrix)
library(psych)

n <- 1000
P <- 100
NUM_SIM <- 50

block <- 10
mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
Sigma <-   bdiag(mat, mat, mat, mat, mat, mat, mat, mat, mat, mat)
corSigma <- cov2cor(Sigma)
icorSigma <- cov2cor(solve(corSigma))

frob_vals <- matrix(0, NUM_SIM, 7)

for(m in 1:NUM_SIM){

  data <- MASS::mvrnorm(n,rep(0,P),Sigma)
  S <- cov(data, method = "pearson")

  ###  GLASSO

  nr  <- 100
  rho <- c(0.01, 0.1, 0.5, 1)
  a_list <- list()
  for(i in 1:length(rho)){
    a_list[[i]] <- glasso::glasso(S,rho[i])
  }

  ## Strimmer Shafer

  strimmer_sample <- cov2cor(corpcor::invcor.shrink(data))

  ## CorShrink

  cov_sample_ML <-  CorShrinkData(data, sd_boot = FALSE,
                                  ash.control = list())
  icor_shrink  <- cov2cor(solve(cov_sample_ML$cor))


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

save(frob_vals, file = paste0("hub_cmd_boot_n_", n, "_P_", P, "_results_icor.rda"))

