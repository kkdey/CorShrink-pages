
##############          eigenvalues distribution (Hub matrix)          ##################

library(ggplot2)
library(corpcor)
library(CorShrink)
library(gridExtra)
library(glasso)
library(Matrix)
library(psych)

gg <- list()

num_samp <- c(10, 50, 100, 1000)
block <- 10
mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
Sigma <-   bdiag(mat, mat, mat, mat, mat, mat, mat, mat, mat, mat)
corSigma <- cov2cor(Sigma)

saved_result <- list()

for(num_iter in 1:length(num_samp)){
  n <- num_samp[num_iter]
  P <- 100
  tmp <- c()

  for(num_sim in 1:30){
    data <- MASS::mvrnorm(n,rep(0,P),corSigma)
    S <- cov(data)


    ###  GLASSO

    nr  <- 100
    rho <- c(1e-04, 0.01, 1)
    a_list <- list()
    for(i in 1:length(rho)){
      a_list[[i]] <- glasso::glasso(S,rho[i])
    }

    ### Strimmer Shafer

    strimmer_sample <- corpcor::cor.shrink(data)


    ### CorShrink

    cov_sample_ML <-  CorShrinkData(data, sd_boot = FALSE,
                                    image = "null",
                                    ash.control = list())

    pdsoft_sample <- PDSCE::pdsoft.cv(data, tolin = 1e-04, tolout = 1e-04)
    pdcor <- cov2cor(pdsoft_sample$sigma)

    num <- 20
    tmp <- rbind(tmp, sqrt(c(eigen(S)$values[1:num], eigen(strimmer_sample)$values[1:num],
                         eigen(cov_sample_ML$cor)$values[1:num],
                         eigen(pdcor)$values[1:num],
                         eigen(cov2cor(a_list[[1]]$w))$values[1:num],
                         eigen(cov2cor(a_list[[2]]$w))$values[1:num],
                         eigen(cov2cor(a_list[[3]]$w))$values[1:num],
                         eigen(corSigma)$values[1:num])))
    cat("We are at iteration", num_sim, "\n")
  }

  tmp_mean <- apply(tmp, 2, mean)
  tmp_sd <- apply(tmp, 2, sd)

  saved_result[[num_iter]] <- list("mean" = tmp_mean, "sd" =  tmp_sd)
  cat("We are at iteration", num_iter, "\n")

}


save(saved_result, file = "hub_sqrt_eigenvalues_distribution.rda")

