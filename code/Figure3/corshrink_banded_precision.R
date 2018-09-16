

##############   CorShrink on the banded precision matrices  ###################


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

frob_vals <- matrix(0, NUM_SIM, 10)

for(m in 1:NUM_SIM){

  data <- MASS::mvrnorm(n,rep(0,P),corSigma)
  S <- cov(data)

  ###  GLASSO

  nr  <- 100
  rho <- c(1e-02, 0.01, 1, 5, 10, 50)
  a_list <- list()
  for(i in 1:length(rho)){
    a_list[[i]] <- glasso::glasso(S,rho[i])
  }
  ## Strimmer Shafer

  strimmer_sample <- corpcor::cor.shrink(data)

  ## CorShrink

  cov_sample_ML <-  CorShrinkData(data, sd_boot = FALSE, image = "null",
                                  ash.control = list())

  pdsoft_sample <- PDSCE::pdsoft.cv(data, tolin = 1e-04, tolout = 1e-04)

  # frob_S <- mean(as.matrix(cov2cor(S)) - as.matrix(cov2cor(corSigma)))^2
  frob_S <- 1 - (tr(as.matrix(cov2cor(S)%*%cov2cor(corSigma))))/(norm(cov2cor(S), type = "F")* norm(corSigma, type = "F"))

  # frob_glasso <- mean(as.matrix(cov2cor(a$w)) -as.matrix(cov2cor(corSigma)))^2

  # frob_strimmer <- mean(as.matrix(cov2cor(strimmer_sample)) - as.matrix(cov2cor(corSigma)))^2
  frob_strimmer <- 1 - (tr(as.matrix(cov2cor(strimmer_sample[1:P,1:P])%*%cov2cor(corSigma))))/(norm(cov2cor(strimmer_sample[1:P,1:P]), type = "F")* norm(corSigma, type = "F"))

  # frob_corshrink <- mean(as.matrix(cov2cor(cov_sample_ML$ash_cor_PD)) - as.matrix(cov2cor(corSigma)))^2
  frob_corshrink <- 1 - (tr(as.matrix(cov2cor(cov_sample_ML$cor)%*%cov2cor(corSigma))))/(norm(cov2cor(cov_sample_ML$cor), type = "F")* norm(corSigma, type = "F"))

  frob_pdsoft <- 1 - (tr(as.matrix(cov2cor(pdsoft_sample$sigma)%*%cov2cor(corSigma))))/(norm(cov2cor(pdsoft_sample$sigma), type = "F")* norm(corSigma, type = "F"))

  frob_glasso_1 <- 1 - (tr(as.matrix(cov2cor(a_list[[1]]$w)%*%cov2cor(corSigma))))/(norm(cov2cor(a_list[[1]]$w), type = "F")* norm(corSigma, type = "F"))
  frob_glasso_2 <- 1 - (tr(as.matrix(cov2cor(a_list[[2]]$w)%*%cov2cor(corSigma))))/(norm(cov2cor(a_list[[2]]$w), type = "F")* norm(corSigma, type = "F"))
  frob_glasso_3 <- 1 - (tr(as.matrix(cov2cor(a_list[[3]]$w)%*%cov2cor(corSigma))))/(norm(cov2cor(a_list[[3]]$w), type = "F")* norm(corSigma, type = "F"))
  frob_glasso_4 <- 1 - (tr(as.matrix(cov2cor(a_list[[4]]$w)%*%cov2cor(corSigma))))/(norm(cov2cor(a_list[[4]]$w), type = "F")* norm(corSigma, type = "F"))
  frob_glasso_5 <- 1 - (tr(as.matrix(cov2cor(a_list[[5]]$w)%*%cov2cor(corSigma))))/(norm(cov2cor(a_list[[5]]$w), type = "F")* norm(corSigma, type = "F"))
  frob_glasso_6 <- 1 - (tr(as.matrix(cov2cor(a_list[[6]]$w)%*%cov2cor(corSigma))))/(norm(cov2cor(a_list[[5]]$w), type = "F")* norm(corSigma, type = "F"))

  frob_vals[m, ] <- c(frob_S, frob_strimmer, frob_corshrink, frob_pdsoft,
                      frob_glasso_1, frob_glasso_2, frob_glasso_3, frob_glasso_4,
                      frob_glasso_5, frob_glasso_6)
  cat("We are at simulation", m, "\n")

}

save(frob_vals, file = paste0("banded_precision_n_", n, "_P_", P, "_results.rda"))
