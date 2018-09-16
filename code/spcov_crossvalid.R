

###############   spcov cross validation   ##########################

library(spcov)

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

Pat <- matrix(1, P, P)
diag(Pat) <- 0

lam <- 0.06
step.size <- 100
tol <- 1e-4

cross_set <- 4

lambda_set <- c (0.01, 0.02, 0.05, 0.08, 0.1, 0.2, 0.5, 0.7, 1, 2)
arr_lambda <- array(0, 10)
for(num_lambda in 1:10){
  arr <- array(0, 10)
  for(NUM_K in 1:10){
    cross_train <- sample(1:n, 4, replace = FALSE)
    cross_test <- setdiff(1:n, cross_train)
    cormat_train <- cor(data[cross_train,], use = "pairwise.complete.obs")
    cormat_test <- cor(data[cross_test,], use = "pairwise.complete.obs")
    mm <- spcov(Sigma=cormat_train+0.1*diag(1,P), S=cormat_train+0.1*diag(1,P), lambda=lambda_set[num_lambda] * Pat,
                step.size=step.size, n.inner.steps=200, thr.inner=0, tol.outer=tol, trace=1)
    arr[NUM_K] = -log(det(mm$Sigma)) - psych::tr(cormat_test %*% solve(mm$Sigma))
    cat("We are at iteration", NUM_K, "\n")
  }

  arr_lambda[num_lambda] = mean(arr)
  cat("We are at iteration", num_lambda, "\n")
}
