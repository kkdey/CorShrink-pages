

#########  corshrink partial correlation hub  ###################


library(CorShrink)
library(huge)
library(corpcor)
library(ggcorrplot)
library(Matrix)
library(psych)

n <- 10
P <- 100
NUM_SIM <- 30

block <- 10
mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
Sigma <-   bdiag(mat, mat, mat, mat, mat, mat, mat, mat, mat, mat)
corSigma <- pcor2cor(as.matrix(Sigma))
pcorSigma <- cor2pcor(corSigma)

ggcorrplot(as.matrix(corSigma))
ggcorrplot(pcorSigma)

frob_vals <- matrix(0, NUM_SIM, 4)

for(m in 1:NUM_SIM){

  data <- MASS::mvrnorm(n,rep(0,P),Sigma)
  S <- cor(data)
  pS <- cor2pcor(S)

  ###  GLASSO

  nr  <- 100
  rho <- seq(0.1,10,length=nr)
  bic <- rho
  S <- cov(data)
  for(j in 1:nr){
    a       <- glasso::glasso(S,rho[j])
    p_off_d <- sum(a$wi!=0 & col(S)<row(S))
    bic[j]  <- -2*(a$loglik) + p_off_d*log(n)
  }
  best <- which.min(bic)

  a <- glasso::glasso(S,rho[best])
  pglasso <- cor2pcor(cov2cor(a$w))

  ## Strimmer Shafer

  strimmer_sample <- corpcor::pcor.shrink(data)

  ## CorShrink

  cov_sample_ML <-  CorShrinkData(data, sd_boot = FALSE, nboot = 300, type = "pcor",
                                  image = "both", ash.control = list())


  # frob_S <- mean(as.matrix(cov2cor(S)) - as.matrix(cov2cor(corSigma)))^2
  frob_S <- 1 - (tr(as.matrix(pS%*%cov2cor(pcorSigma))))/(norm(pS, type = "F")* norm(pcorSigma, type = "F"))

  # frob_glasso <- mean(as.matrix(cov2cor(a$w)) -as.matrix(cov2cor(corSigma)))^2
  frob_glasso <- 1 - (tr(as.matrix(pglasso%*%pcorSigma)))/(norm(pglasso, type = "F")* norm(pcorSigma, type = "F"))

  # frob_strimmer <- mean(as.matrix(cov2cor(strimmer_sample)) - as.matrix(cov2cor(corSigma)))^2
  frob_strimmer <- 1 - (tr(as.matrix(strimmer_sample[1:P,1:P]%*%pcorSigma)))/(norm(strimmer_sample[1:P,1:P], type = "F")* norm(pcorSigma, type = "F"))

  # frob_corshrink <- mean(as.matrix(cov2cor(cov_sample_ML$ash_cor_PD)) - as.matrix(cov2cor(corSigma)))^2
  frob_corshrink <- 1 - (tr(as.matrix(cov2cor(cov_sample_ML$ash_cor_PD)%*%cov2cor(pcorSigma))))/(norm(cov2cor(cov_sample_ML$ash_cor_PD), type = "F")* norm(pcorSigma, type = "F"))
 # frob_corshrink <- 1 - (tr(as.matrix(cov2cor(cov_sample_ML$ash_cor_PD)%*%cov2cor(corSigma))))/(norm(cov2cor(cov_sample_ML$ash_cor_PD), type = "F")* norm(corSigma, type = "F"))


  frob_vals[m, ] <- c(frob_S, frob_glasso, frob_strimmer, frob_corshrink)
  cat("We are at simulation", m, "\n")

}

save(frob_vals, file = "../output/hub/hub_cmd_partial_boot_n_10_P_100_results.rda")

df <- data.frame("method" = c(rep("empirical", dim(frob_vals)[1]),
                              rep("glasso", dim(frob_vals)[1]),
                              rep("corpcor", dim(frob_vals)[1]),
                              rep("corshrink", dim(frob_vals)[1])),
                 "distance" = log(c(frob_vals[,1], frob_vals[,2],
                                    frob_vals[,3], frob_vals[,4])))

p <- ggplot(df, aes(method, distance)) + xlab("") + ylab("log(distance)") + ggtitle("n=10, p = 100")
p1 <- p + geom_boxplot()
print(p1)
