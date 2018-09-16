

###############   adaptive shrinkage asymptotics  ############################

library(ashr)
beta <- rnorm(100, 0, 1)
betahat <- sapply(beta, function(x) return(rnorm(1, x, 0.01)))
ash_betahat <- ash(betahat, rep(0.01, length(betahat)))
mean((betahat - beta)^2)
mean((ash_betahat$result$PosteriorMean - beta)^2)
