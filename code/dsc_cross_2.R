
set.seed(100)

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

n <- 50
P <- 100
ll <- DM_toeplitz(n=n, P=P)
data <- rbind(ll$Xtrain, ll$Xtest)
Sigma <- ll$Sigma
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
glasso_sample <- glasso::glasso(covmat, rho = rho_array[which.max(score)])
cov_sample_ML <-  CorShrinkMatrix(cov2cor(covmat), nsamp = matrix(n, P, P),
                                  ash.control = list(mixcompdist = "normal",
                                                     nullweight = 10))

lam <- 0.2
step.size <- 10
tol <- 0.01
covmat <- cov(data)

lambda_set <- c (0.02, 0.05, 0.08, 0.1, 0.2, 0.5, 0.7, 1 )

score_spcor <- 1/abs(covmat)
diag(score_spcor) = 0
score_spcor[score_spcor > 10^5] = 10^5
score_spcor[score_spcor < -10^5] = -10^5

arr_lambda <- array(0, length(lambda_set))
for(num_lambda in 1:length(lambda_set)){
  arr <- array(0, 1)
  for(NUM_K in 1:1){
    cross_train <- sample(1:n, round(0.2*n), replace = FALSE)
    cross_test <- setdiff(1:n, cross_train)
    cormat_train <- cor(data[cross_train,], use = "pairwise.complete.obs")
    cormat_test <- cor(data[cross_test,], use = "pairwise.complete.obs")
    mm <- spcov(Sigma=cormat_train+0.1*diag(1,P),
                S=cormat_train+0.1*diag(1,P),
                lambda=lambda_set[num_lambda] * score_spcor,
                step.size=step.size, n.inner.steps=20,
                thr.inner=0, tol.outer=0.01, trace=1)
    arr[NUM_K] = -log(det(mm$Sigma)) -
                      psych::tr(cormat_test %*% solve(mm$Sigma))
  }

  arr_lambda[num_lambda] = mean(arr)
  cat("We are at iteration", num_lambda, "\n")
}

mm <- spcov(Sigma=covmat+0.1*diag(1,P), S=covmat+0.1*diag(1,P),
             lambda=lambda_set[which.max(arr_lambda)] * score_spcor,
            step.size=step.size, n.inner.steps=20,
            thr.inner=0, tol.outer=tol, trace=1)
spcor_mat <- cov2cor(mm$Sigma)


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


image(as.matrix(spcor_mat),
      col=col, main=paste0("spcov: "), cex.main=2,
      xaxt = "n", yaxt = "n", zlim=c(-1,1))



