

################  CorShrink banded precision (pcor)  #############################


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
pcorSigma <- cor2pcor(corSigma)


frob_vals <- matrix(0, NUM_SIM, 8)

frob_vals <- matrix(0, NUM_SIM, 8)

for(m in 1:NUM_SIM){

  data <- MASS::mvrnorm(n,rep(0,P),Sigma)
  S <- cov(data, method = "pearson")

  ###  GLASSO

  nr  <- 100
  rho <- c(0.01, 0.1, 0.5, 1)
  a_list <- list()

  for(i in 1:length(rho)){
    tmp <- glasso::glasso(S,rho[i])
    wi <- tmp$wi
    pcor_from_wi <- -cov2cor(wi)
    diag(pcor_from_wi) <- rep(1, dim(pcor_from_wi)[1])
    a_list[[i]] <- pcor_from_wi
  }

  ## Strimmer Shafer

  strimmer_sample <- corpcor::pcor.shrink(data)

  ## CorShrink

  cov_sample_ML <-  pCorShrinkData(data, nboot = 50)
  pcor_shrink  <- cov_sample_ML$cor

  pcorS <- cor2pcor(cov2cor(S))

  pdsoft_sample <- PDSCE::pdsoft.cv(data, tolin = 1e-04, tolout = 1e-04)
  pcor_pdsoft_sample <- cor2pcor(cov2cor(pdsoft_sample$sigma))

  frob_S <-   mean(sqrt((as.matrix(pcorS) - as.matrix(pcorSigma))^2))
  frob_strimmer <-   mean(sqrt((as.matrix(strimmer_sample[1:P,1:P]) - as.matrix(pcorSigma))^2))
  frob_corshrink <- mean(sqrt((pcor_shrink - as.matrix(pcorSigma))^2))
  frob_glasso_1 <- mean(sqrt((as.matrix(a_list[[1]]) - as.matrix(pcorSigma))^2))
  frob_glasso_2 <- mean(sqrt((as.matrix(a_list[[2]]) - as.matrix(pcorSigma))^2))
  frob_glasso_3 <- mean(sqrt((as.matrix(a_list[[3]]) - as.matrix(pcorSigma))^2))
  frob_glasso_4 <- mean(sqrt((as.matrix(a_list[[4]]) - as.matrix(pcorSigma))^2))
  frob_pdsoft <- mean(sqrt((as.matrix(pcor_pdsoft_sample) - as.matrix(pcorSigma))^2))


  frob_vals[m, ] <- c(frob_S, frob_strimmer, frob_corshrink,  frob_pdsoft,
                      frob_glasso_1,
                      frob_glasso_2,  frob_glasso_3,  frob_glasso_4)
  cat("We are at simulation", m, "\n")

}

save(frob_vals, file = paste0("banded_precision_n_", n, "_P_", P, "_results_pcor.rda"))


