

###########  Combine GLASSO + CorShrink  #####################

library(CorShrink)
library(glasso)
library(Matrix)
library(pracma)
library(corrplot)
library(corpcor)

#################  Choice of n and P  ############################
n <- 1000
P <- 100


#################  Generate the hub correlation matrix  #########################
block <- 10
mat <- 0.3*diag(1,block) + 0.7*rep(1,block) %*% t(rep(1, block))
Sigma <-   bdiag(mat, mat, mat, mat, mat, mat, mat, mat, mat, mat)
corSigma <- cov2cor(Sigma)
pcorSigma <- cor2pcor(corSigma)

##################  Which samples are used for GLASSO fit and which for CorShrink

glasso_indices <- sample(1:n, floor(n/2), replace = FALSE)
corshrink_indices <- setdiff(1:n, glasso_indices)

#################   Data generation  ##############################

data <- MASS::mvrnorm(n,rep(0,P),Sigma)
S <- cov(data, method = "pearson")


#################   fit GLASSO on partial data  #######################

S_glasso <- cov(data[glasso_indices, ], method = "pearson")
out <- glasso(S_glasso, 0.1) ### lambda = 0.1
L <- sqrtm(out$wi)$B
data_new <- t(L %*% t(data[corshrink_indices, ]))
corshrink_out <- CorShrinkData(data_new)$cor
inv_est <- L %*% corshrink_out %*% t(L)
pcor_est <- -cov2cor(inv_est)
diag(pcor_est) <- rep(1, dim(pcor_est)[1])

col2 <- c("blue", "white", "red")
corrplot(pcor_est, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

#############  Applyting GLASSO directly  ##############

glasso1 <- glasso(S, 0.1)
pcor_est2 <- -cov2cor(glasso1$wi)
diag(pcor_est2) <- rep(1, dim(pcor_est2)[1])

col2 <- c("blue", "white", "red")
corrplot(pcor_est2, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.2, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")

############  Frobenius distance  ###################

mean(sqrt((as.matrix(pcor_est) - as.matrix(pcorSigma))^2))
mean(sqrt((as.matrix(pcor_est2) - as.matrix(pcorSigma))^2))
