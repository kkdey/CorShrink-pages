

###################   DSC  cross validation examples   ############################

set.seed(100)

DM_diagonal = function(n,P){
  library("MASS")
  index1=sort(sample(seq(1:n),(n/2)))
  index2=seq(1:n)[-index1]
  Sigma = diag(rchisq(P,3))
  data = mvrnorm(n,rep(0,P),Sigma)
  Xtest = data[index2,]
  Xtrain = data[index1,]
  Omega = diag((1/rchisq(P,3)))
  return(list(Xtrain = Xtrain, Xtest = Xtest, Omega = Omega))
}


n <- 10
P <- 100
ll <- DM_diagonal(n=n, P=P)
data <- rbind(ll$Xtrain, ll$Xtest)
Sigma <- solve(ll$Omega)
corSigma <- cov2cor(Sigma)
K <- 20
rho_array <- c(0.005, 0.05, 0.1, 0.2, 0.5, 0.8, 1, 2, 5, 10, 100)
covmat <- cov(data)
score<- array(0, length(rho_array))
for(num in 1:length(rho_array)){
  temp <- array(0, K)
  for(k in 1:K){
    cross_train <- sample(1:n, round(0.2*n), replace = FALSE)
    cross_test <- setdiff(1:n, cross_train)
    covmat_train <- cov(data[cross_train,], use = "pairwise.complete.obs")
    covmat_test <- cov(data[cross_test,], use = "pairwise.complete.obs")
    glasso_train <- as.matrix(cov2cor(glasso::glasso(covmat_train, rho = rho_array[num])$w))
    temp[k] <- -log(det(glasso_train)) - psych::tr(covmat_test %*% solve(glasso_train))
  }
  cat("We are at num", num, "\n")
  score[num] <- mean(temp)
}

strimmer_sample <- corpcor::cov.shrink(data)
glasso_sample <- glasso::glasso(cov_mat, rho = 100)
cov_sample_ML <-  CorShrinkMatrix(cov2cor(cov_mat), nsamp = matrix(n, P, P),
                                  ash.control = list(mixcompdist = "normal", nullweight = 100))

col=c(rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)))

par(mfrow = c(3,2))

image(as.matrix(corSigma),
      col=col, main=paste0("pop corr: "), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))

image(as.matrix(cov2cor(covmat)),
      col=col, main=paste0("sample corr: "), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))


image(as.matrix(cov2cor(strimmer_sample)),
      col=col, main=paste0("shafer strimmer: "), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))

image(as.matrix(cov2cor(glasso_sample$w)),
      col=col, main=paste0("glasso : "), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))


image(as.matrix(cov_sample_ML$ash_cor_PD),
      col=col, main=paste0("corshrink: "), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))





