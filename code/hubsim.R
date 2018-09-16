
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

