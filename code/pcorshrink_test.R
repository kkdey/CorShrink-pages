

###################   test pCorShrink ()   #######################

library(CorShrink)
library(glasso)
library(corpcor)
library(Matrix)
library(psych)
library(corrplot)


n <- 1000
P <- 100
block <- 10
mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
Sigma <-   bdiag(mat, mat, mat, mat, mat, mat, mat, mat, mat, mat)
corSigma <- cov2cor(Sigma)
icorSigma <- solve(corSigma)
pcorSigma <- corpcor::cor2pcor(corSigma)


data <- MASS::mvrnorm(n,rep(0,P),Sigma)
out <- pCorShrinkData(data)

col2 <- c("blue", "white", "red")
corrplot(out$cor, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")


