

######################  Eigenvalues distribution    #############################

library(ggplot2)
library(corpcor)
library(CorShrink)
library(gridExtra)
library(glasso)

gg <- list()

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


num_samp <- c(10, 50, 100, 1000)
saved_result <- list()

for(num_iter in 1:length(num_samp)){

  n <- num_samp[num_iter]
  P <- 100
  tmp <- c()

  for(num_sim in 1:20){
    ll <- DM_toeplitz(n=n, P=P)
    data <- rbind(ll$Xtrain, ll$Xtest)
    Sigma <- ll$Sigma
    corSigma <- cov2cor(Sigma)

    S <- cov2cor(cov(data))

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
                                    ash.control = list())

    pdsoft_sample <- PDSCE::pdsoft.cv(data, tolin = 1e-04, tolout = 1e-04)
    pdcor <- cov2cor(pdsoft_sample$sigma)

    num <- 20

    num <- 20
    tmp <- rbind(tmp, sqrt(c(eigen(S)$values[1:num], eigen(strimmer_sample)$values[1:num],
                             eigen(cov_sample_ML$cor)$values[1:num],
                             eigen(pdcor)$values[1:num],
                             eigen(cov2cor(a_list[[1]]$w))$values[1:num],
                             eigen(cov2cor(a_list[[2]]$w))$values[1:num],
                             eigen(cov2cor(a_list[[3]]$w))$values[1:num],
                             eigen(corSigma)$values[1:num])))
    cat("We are at simulation", num_sim, "\n")
  }

  tmp_mean <- apply(tmp, 2, mean)
  tmp_sd <- apply(tmp, 2, sd)
  saved_result[[num_iter]] <- list("mean"= tmp_mean, "sd" = tmp_sd)
  cat("We are at iteration", num_iter, "\n")

}

save(saved_result, file = "toeplitz_sqrt_eigenvalues_distribution.rda")


# grid.arrange(gg[[1]]$plot, gg[[2]]$plot, gg[[3]]$plot,
#              gg[[4]]$plot, nrow = 4, ncol = 1)


