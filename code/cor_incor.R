library(CorShrink)
library(glasso)
library(corpcor)
library(Matrix)
library(psych)
library(corrplot)

n <- 10
P <- 100

block <- 10
mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
Sigma <-   bdiag(mat, mat, mat, mat, mat, mat, mat, mat, mat, mat)
corSigma <- cov2cor(Sigma)
icorSigma <- solve(corSigma)

col2 <- c("blue", "white", "red")
corrplot((as.matrix(corSigma)), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

pcorSigma <- -cov2cor(as.matrix(icorSigma))
diag(pcorSigma) <- rep(1, dim(pcorSigma)[1])

col2 <- c("blue", "white", "red")
corrplot(pcorSigma, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")


data <- MASS::mvrnorm(n,rep(0,P),Sigma)
S <- cov(data, method = "pearson")
corrplot(cov2cor(S), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")


############  pseudo-inverse of S  ################################

pseudo_S <- pseudoinverse(cov2cor(S), tol = 1e-10)
pcor_S <- -cov2cor(pseudo_S)
diag(pcor_S) <- rep(1, dim(pcor_S)[1])


corrplot(pcor_S, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

########### bagged pseudo-inverse cor of S  ##########################

bags <- 100
pseudo_S_bag <- 0
for(num in 1:bags){
  data_samp <- data[sample(1:dim(data)[1], dim(data)[1], replace = TRUE), ]
  S_bag <- cor(data)
  tmp <- -cov2cor(pseudoinverse(S_bag, tol = 1e-10))
  diag(tmp) <- rep(1, dim(tmp)[1])
  pseudo_S_bag =  pseudo_S_bag  + tmp
}

pS2 <- pseudo_S_bag/bags
corrplot(pS2, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")
out <- cor.fit.mixture(pS2)

cov_sample_ML <-  CorShrinkMatrix(pS2, nsamp = 30,
                                  ash.control = list(mixcompdist = "halfuniform"))
corrplot(cov2cor(cov_sample_ML$cor), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

corrplot(cov2cor(pS2), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")


########### pseudo-inverse bagged cor of S  ##########################

bags <- 100
S_bag <- 0
for(num in 1:bags){
  data_samp <- data[sample(1:dim(data)[1], dim(data)[1], replace = TRUE), ]
  tmp <- cor(data)
  S_bag =  S_bag + tmp
}

S_bag <- S_bag/bags
pseudo_S_bag <- pseudoinverse(S_bag, tol = 1e-10)
pS3 <- -cov2cor(pseudo_S_bag)
diag(pS3) <- rep(1, dim(pS3)[1])
corrplot(pS3, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")
out <- cor.fit.mixture(pS3)







pseudo_S_bag <- pseudo_S_bag/100
diag(pseudo_S_bag) <- rep(1, dim(pseudo_S_bag)[1])
corrplot(pseudo_S_bag, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

cor2pcor(corSigma)
a2 <- -cov2cor(solve(corSigma))
cov_sample_ML <-  CorShrinkMatrix(pseudo_S_bag, nsamp = 12, ash.control = list())
corrplot(cov2cor(cov_sample_ML$cor), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

cov_sample_ML <-  CorShrinkMatrix(pcor_S, nsamp = 3, ash.control = list())
corrplot(cov2cor(cov_sample_ML$cor), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")



cor.fit.mixture <- function(r)
{
  if ( any(r > 1) || any(r < -1) )
    stop("Data out of range: input correlations must be in [-1; 1]")
  zcor <- atanh(r)
  zcorvec <- zcor[lower.tri(zcor)]
  out <- locfdr::locfdr(zcorvec)
  eta0 <- as.double( out$fp0[1,3] ) # p0
  sigma <- as.double( out$fp0[1,2] ) # sig
  kappa <- 1/(sigma*sigma) + 2 # Fisher's rule
  prob.nonzero <- 1-out$fdr
  return( list(kappa=kappa, eta0=eta0, prob.nonzero=prob.nonzero) )
}

out <- cor.fit.mixture(pseudo_S_bag)

zcor <- atanh(pseudo_S_bag)
plot(density(zcor))
corrplot(pseudo_S_bag, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")




df=7
plot.locfdr=0
r <- pcor_S



strimmer_sample <- solve(corpcor::cor.shrink(data))
est <- -cov2cor(strimmer_sample)
diag(est) <- rep(1, dim(strimmer_sample)[1])
corrplot(est, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")


corrplot(corpcor::cor.shrink(data), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")




## CorShrink

cov_sample_ML <-  CorShrinkData(data, sd_boot = FALSE,
                                ash.control = list())
icor_shrink  <- solve(cov_sample_ML$cor)

plot(eigen(cov_sample_ML$cor)$values, ylim = c(0, 10))
plot(eigen(strimmer_sample)$values, ylim = c(0,10))
plot(eigen(corSigma)$values, ylim = c(0, 10))

plot(eigen(cov2cor(S))$values, ylim = c(0, 10))


plot(eigen(icor_shrink)$values, ylim = c(0, 20))
