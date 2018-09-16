
#############  CorShrink on hub matrix  ###################

library(CorShrink)
library(glasso)
library(corpcor)
library(Matrix)
library(psych)
library(PDSCE)

n <- 50
P <- 100
NUM_SIM <- 50

block <- 10
mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
Sigma <-   bdiag(mat, mat, mat, mat, mat, mat, mat, mat, mat, mat)
corSigma <- cov2cor(Sigma)

frob_vals <- matrix(0, NUM_SIM, 10)

for(m in 1:NUM_SIM){

  data <- MASS::mvrnorm(n,rep(0,P),Sigma)
  S <- cov(data, method = "pearson")

  ###  GLASSO

  nr  <- 100
  rho <- c(1e-04, 0.01, 1, 5, 10, 50)
  a_list <- list()
  for(i in 1:length(rho)){
    a_list[[i]] <- glasso::glasso(S,rho[i])
  }

  ## Strimmer Shafer

  strimmer_sample <- corpcor::cor.shrink(data)

  ## CorShrink

  cov_sample_ML <-  CorShrinkData(data, sd_boot = FALSE,
                                  ash.control = list())

  pdsoft_sample <- PDSCE::pdsoft.cv(data, tolin = 1e-04, tolout = 1e-04)


  frob_S <- mean(sqrt((as.matrix(cov2cor(S)) - as.matrix(cov2cor(corSigma)))^2))

  frob_strimmer <- mean(sqrt((as.matrix(cov2cor(strimmer_sample[1:P,1:P])) - as.matrix(cov2cor(corSigma)))^2))

  frob_corshrink <- mean(sqrt((as.matrix(cov2cor(cov_sample_ML$cor)) - as.matrix(cov2cor(corSigma)))^2))

  frob_pdsoft <- mean(sqrt((as.matrix(cov2cor(pdsoft_sample$sigma)) - as.matrix(cov2cor(corSigma)))^2))

  frob_glasso_1 <- mean(sqrt((as.matrix(cov2cor(a_list[[1]]$w)) - as.matrix(cov2cor(corSigma)))^2))
  frob_glasso_2 <- mean(sqrt((as.matrix(cov2cor(a_list[[2]]$w)) - as.matrix(cov2cor(corSigma)))^2))
  frob_glasso_3 <- mean(sqrt((as.matrix(cov2cor(a_list[[3]]$w)) - as.matrix(cov2cor(corSigma)))^2))
  frob_glasso_4 <- mean(sqrt((as.matrix(cov2cor(a_list[[4]]$w)) - as.matrix(cov2cor(corSigma)))^2))
  frob_glasso_5 <- mean(sqrt((as.matrix(cov2cor(a_list[[5]]$w)) - as.matrix(cov2cor(corSigma)))^2))
  frob_glasso_6 <- mean(sqrt((as.matrix(cov2cor(a_list[[6]]$w)) - as.matrix(cov2cor(corSigma)))^2))



  frob_vals[m, ] <- c(frob_S, frob_strimmer, frob_corshrink, frob_pdsoft,
                      frob_glasso_1, frob_glasso_2, frob_glasso_3, frob_glasso_4,
                      frob_glasso_5, frob_glasso_6)

  cat("We are at simulation", m, "\n")

}

save(frob_vals, file = paste0("hub_frobenius_boot_n_", n, "_P_", P, "_results.rda"))

