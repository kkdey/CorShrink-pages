

###########  corshrink on Toeplitz Pearson, Spearman and Kendall  ######################

library(CorShrink)
library(huge)
library(corpcor)

n <- 10
P <- 100
NUM_SIM <- 10


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


frob_vals <- matrix(0, NUM_SIM, 6)

for(m in 1:NUM_SIM){
  ll <- DM_toeplitz(n=n, P=P)
  data <- rbind(ll$Xtrain, ll$Xtest)
  Sigma <- ll$Sigma
  corSigma <- cov2cor(Sigma)

  ###  GLASSO


  S <- cov(data)

  kendall_noboot <-  CorShrinkData(data, sd_boot = FALSE,
                                  nboot = 200,
                                  cor_method = "kendall")

  kendall_boot <-  CorShrinkData(data, sd_boot = TRUE,
                                   nboot = 200,
                                   cor_method = "kendall")


  spearman_noboot <-  CorShrinkData(data, sd_boot = FALSE,
                                   nboot = 200,
                                   cor_method = "spearman")

  spearman_boot <-  CorShrinkData(data, sd_boot = TRUE,
                                 nboot = 200,
                                 cor_method = "spearman")

  pearson_noboot <-  CorShrinkData(data, sd_boot = FALSE,
                                    nboot = 200,
                                    cor_method = "pearson")

  pearson_boot <-  CorShrinkData(data, sd_boot = TRUE,
                                 nboot = 200,
                                 cor_method = "pearson")


  frob_pearson_boot <- 1 - (tr(as.matrix(cov2cor(pearson_boot$ash_cor_PD)%*%cov2cor(corSigma))))/(norm(cov2cor(pearson_boot$ash_cor_PD), type = "F")* norm(corSigma, type = "F"))
  frob_pearson_noboot <- 1 - (tr(as.matrix(cov2cor(pearson_noboot$ash_cor_PD)%*%cov2cor(corSigma))))/(norm(cov2cor(pearson_noboot$ash_cor_PD), type = "F")* norm(corSigma, type = "F"))
  frob_spearman_boot <- 1 - (tr(as.matrix(cov2cor(spearman_boot$ash_cor_PD)%*%cov2cor(corSigma))))/(norm(cov2cor(spearman_boot$ash_cor_PD), type = "F")* norm(corSigma, type = "F"))
  frob_spearman_noboot <- 1 - (tr(as.matrix(cov2cor(spearman_noboot$ash_cor_PD)%*%cov2cor(corSigma))))/(norm(cov2cor(spearman_noboot$ash_cor_PD), type = "F")* norm(corSigma, type = "F"))
  frob_kendall_boot <- 1 - (tr(as.matrix(cov2cor(kendall_boot$ash_cor_PD)%*%cov2cor(corSigma))))/(norm(cov2cor(kendall_boot$ash_cor_PD), type = "F")* norm(corSigma, type = "F"))
  frob_kendall_noboot <- 1 - (tr(as.matrix(cov2cor(kendall_noboot$ash_cor_PD)%*%cov2cor(corSigma))))/(norm(cov2cor(kendall_noboot$ash_cor_PD), type = "F")* norm(corSigma, type = "F"))


  frob_vals[m, ] <- c(frob_pearson_boot, frob_pearson_noboot, frob_spearman_boot,
                      frob_spearman_noboot, frob_kendall_boot, frob_kendall_noboot)
  cat("We are at simulation", m, "\n")
}

save(frob_vals, file = "../output/toeplitz/psk_toeplitz_cmd_n_10_P_100_results.rda")
