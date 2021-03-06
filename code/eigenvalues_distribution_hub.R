

##############          eigenvalues distribution (Hub matrix)          ##################

library(ggplot2)
library(corpcor)
library(CorShrink)
library(gridExtra)

gg <- list()

num_samp <- c(10, 50, 100, 1000)
block <- 10
mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
Sigma <-   bdiag(mat, mat, mat, mat, mat, mat, mat, mat, mat, mat)
corSigma <- cov2cor(Sigma)

for(num_iter in 1:length(num_samp)){
  n <- num_samp[num_iter]
  P <- 100
  tmp <- c()

  for(num_sim in 1:20){
    data <- MASS::mvrnorm(n,rep(0,P),corSigma)
    S <- cov(data)


    ###  GLASSO

    nr  <- 100
    rho <- seq(0.1,10,length=nr)
    bic <- rho
    for(j in 1:nr){
      a       <- glasso::glasso(S,rho[j])
      p_off_d <- sum(a$wi!=0 & col(S)<row(S))
      bic[j]  <- -2*(a$loglik) + p_off_d*log(n)
    }
    best <- which.min(bic)

    a <- glasso::glasso(S,rho[best])
    glasso_cor <- cov2cor(a$w)

    ### Strimmer Shafer

    strimmer_sample <- corpcor::cor.shrink(data)


    ### CorShrink

    cov_sample_ML <-  CorShrinkData(data, sd_boot = FALSE, optmethod = "mixEM",
                                    image.control = list(x.cex = 0.3, y.cex = 0.3),
                                    ash.control = list())

    num <- 20
    tmp <- rbind(tmp, (c(eigen(S)$values[1:num], eigen(strimmer_sample)$values[1:num],
                         eigen(glasso_cor)$values[1:num], eigen(cov_sample_ML$ash_cor_only)$values[1:num],
                         eigen(corSigma)$values[1:num])))
    cat("We are at iteration", num_sim, "\n")
  }

  tmp_mean <- apply(tmp, 2, mean)
  tmp_sd <- apply(tmp, 2, sd)

  saved_result[[num_iter]] <- list("mean" = tmp_mean, "sd" =  tmp_sd)
  cat("We are at iteration", num_iter, "\n")

}


save(saved_result, file = "hub_eigenvalues_distribution.rda")

