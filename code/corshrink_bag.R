

###############  CorShrink  + Bagging partial correlation  ############################

library(CorShrink)
library(glasso)
library(corpcor)
library(Matrix)
library(psych)
library(corrplot)

cor2pcor1 <- function(data){
  return(cor2pcor(cor(data)))
}

cor2pcor2 <- function(data, number_bags = 100){
  pS_bag <- 0
  for(num in 1:number_bags){
    data_samp <- data[sample(1:dim(data)[1], dim(data)[1], replace = TRUE), ]
    S_bag <- cor(data_samp)
    pS_bag <- pS_bag + corpcor::cor2pcor(S_bag)
  }
  pS_bag <- pS_bag/number_bags
  return(pS_bag)
}

cor2pcor3 <- function(data, number_bags = 100){
  S_bag <- 0
  for(num in 1:number_bags){
    data_samp <- data[sample(1:dim(data)[1], dim(data)[1], replace = TRUE), ]
    S_bag <- S_bag + cor(data_samp)
  }
  pS_bag <- cor2pcor((S_bag/number_bags))
  return(pS_bag)
}



n <- 10
P <- 100



block <- 10
mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
Sigma <-   bdiag(mat, mat, mat, mat, mat, mat, mat, mat, mat, mat)
corSigma <- cov2cor(Sigma)
icorSigma <- solve(corSigma)
pcorSigma <- corpcor::cor2pcor(corSigma)

col2 <- c("blue", "white", "red")
corrplot(cov2cor(as.matrix(icorSigma)), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")


data <- MASS::mvrnorm(n,rep(0,P),Sigma)




S <- cov(data, method = "pearson")
pS <- corpcor::cor2pcor(cov2cor(S))

S1 <- cov2cor(cov(data))
corrplot(cov2cor(pseudoinverse(S1)), diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")




out <- pseudoinverse(S1, tol = 0.01)



bags <- 30
pS_bag_summ <- array(0, c(dim(pS)[1], dim(pS)[2], bags))
for(num in 1:bags){
  data_samp <- data[sample(1:dim(data)[1], dim(data)[1], replace = TRUE), ]
  pS_bag <- cor2pcor1(data_samp)
  pS_bag_summ[,,num] <- pS_bag
}

tmp <- atanh(pS_bag_summ)
tmp[tmp==Inf] = NA
sd_pS_bag <- apply(tmp, c(1,2), function(x) return(sd(x)))

out <- CorShrinkMatrix(pS, zscore_sd = sd_pS_bag)

corrplot(pS, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

corrplot(pcorSigma, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")



corrplot(out$cor, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")



strimmer_sample <- corpcor::pcor.shrink(data)

corrplot(strimmer_sample, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

corrplot(pS, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")


mean(sqrt((as.matrix(out$cor) - as.matrix(pcorSigma))^2))
mean(sqrt((as.matrix(pS) - as.matrix(pcorSigma))^2))
mean(sqrt((as.matrix(strimmer_sample) - as.matrix(pcorSigma))^2))
mean(sqrt((0 - as.matrix(pcorSigma))^2))
