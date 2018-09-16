

################  CorShrink banded precision (pcor)  #############################


library(CorShrink)
library(corpcor)
library(glasso)
library(Matrix)
library(psych)

n <- 10
P <- 100
NUM_SIM <- 50

diags <- list()
diags[[1]] <- rep(1, 100)
diags[[2]] <- rep(-0.5, 100)
Kinv <- bandSparse(100, k = -(0:1), diag = diags, symm = TRUE)
K <- solve(Kinv)
corSigma <- cov2cor(K)
pcorSigma <- cor2pcor(corSigma)


col=c(rep(rgb(0,1,0), 10000),rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)), rep(rgb(1,0,0),10000))

image(as.matrix(corSigma),
      col=col, main=paste0("pop corr matrix"), cex.main=1,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))

col=c(rep(rgb(0,1,0), 10000),rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)), rep(rgb(1,0,0),10000))

image(as.matrix(Kinv),
      col=col, main=paste0("pop inv corr matrix"), cex.main=1,
      xaxt = "n", yaxt = "n", zlim = c(-1,1))



frob_vals <- matrix(0, NUM_SIM, 4)

for(m in 1:NUM_SIM){

  data <- MASS::mvrnorm(n,rep(0,P),corSigma)
  S <- cov(data)

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

  ## Strimmer Shafer

  strimmer_sample <- corpcor::pcor.shrink(data)

  ## CorShrink

  cov_sample_ML <-  CorShrinkData(data, sd_boot = TRUE,nboot = 100,
                                  type = "pcor", cor_method = "pearson",
                                  image = "null",
                                  ash.control = list())
  pcor_shrink  <- cov_sample_ML$ash_cor_PD



  # frob_S <- mean(as.matrix(cov2cor(S)) - as.matrix(cov2cor(corSigma)))^2
  frob_S <- 1 - (tr(as.matrix(cor2pcor(cov2cor(S))%*%pcorSigma)))/(norm(cor2pcor(cov2cor(S)), type = "F")* norm(pcorSigma, type = "F"))

  # frob_glasso <- mean(as.matrix(cov2cor(a$w)) -as.matrix(cov2cor(corSigma)))^2
  frob_glasso <- 1 - (tr(as.matrix(cor2pcor(cov2cor(a$w))%*%pcorSigma)))/(norm(cor2pcor(cov2cor(a$w)), type = "F")* norm(pcorSigma, type = "F"))

  # frob_strimmer <- mean(as.matrix(cov2cor(strimmer_sample)) - as.matrix(cov2cor(corSigma)))^2
  frob_strimmer <- 1 - (tr(as.matrix(strimmer_sample[1:P,1:P]%*%pcorSigma)))/(norm(strimmer_sample[1:P,1:P], type = "F")* norm(pcorSigma, type = "F"))

  # frob_corshrink <- mean(as.matrix(cov2cor(cov_sample_ML$ash_cor_PD)) - as.matrix(cov2cor(corSigma)))^2
  frob_corshrink <- 1 - (tr(as.matrix(pcor_shrink%*%pcorSigma)))/(norm(pcor_shrink, type = "F")* norm(pcorSigma, type = "F"))

  frob_vals[m, ] <- c(frob_S, frob_glasso, frob_strimmer, frob_corshrink)
  cat("We are at simulation", m, "\n")

}

save(frob_vals, file = paste0("banded_precision_n_", n, "_P_", P, "_results_pcor.rda"))


