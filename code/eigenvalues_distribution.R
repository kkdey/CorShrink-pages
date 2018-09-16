

######################  Eigenvalues distribution    #############################

library(ggplot2)
library(corpcor)
library(CorShrink)
library(gridExtra)

gg <- list()

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


num_samp <- c(10, 50, 100, 1000)

for(num_iter in 1:length(num_samp)){

  n <- num_samp[num_iter]
  P <- 100
  tmp <- c()

  for(num_sim in 1:20){
    ll <- DM_toeplitz(n=n, P=P)
    data <- rbind(ll$Xtrain, ll$Xtest)
    Sigma <- ll$Sigma
    corSigma <- cov2cor(Sigma)

    S <- cov2cor(cov(data))

    ###  GLASSO

    nr  <- 100
    rho <- seq(1e-04,10,length=nr)
    bic <- rho
    for(j in 1:nr){
      a       <- glasso::glasso(S,rho[j])
      p_off_d <- sum(a$wi!=0 & col(S)<row(S))
      bic[j]  <- -2*(a$loglik) + p_off_d*log(n)
    }
    best <- which.min(bic)

    a <- glasso::glasso(S,rho[best])
    glasso_cor <- cov2cor(a$w)

    ### Strimmer Shafer

    strimmer_sample <- corpcor::cor.shrink(data)


    ### CorShrink

    cov_sample_ML <-  CorShrinkData(data, sd_boot = FALSE,
                                    ash.control = list())

    num <- 20

    tmp <- rbind(tmp, (c(eigen(S)$values[1:num], eigen(strimmer_sample)$values[1:num],
                         eigen(glasso_cor)$values[1:num], eigen(cov_sample_ML$ash_cor_only)$values[1:num],
                         eigen(corSigma)$values[1:num])))
    cat("We are at simulation", num_sim, "\n")
  }

  tmp_mean <- apply(tmp, 2, mean)
  # eigs.df <- data.frame ("x" = rep(1:num, 5),
  #                        "y" = tmp_mean,
  #                        "color" = c(rep("empirical", num), rep("corpcor",  num), rep("glasso", num),
  #                                    rep("corshrink", num), rep("true", num)),
  #                        "type" = c(rep("A", 4*num), rep("B", num)))
  #
  # p <- ggplot(eigs.df, aes(x=x, y=y, colour=color, linetype = color)) + geom_line(lty = 1, lwd = 0.7) +
  #   scale_linetype_manual(values = c(rep("solid", 4), rep("dashed", 1))) +
  #   scale_colour_manual(values=c("#000000", "blue", "green", "gold",
  #                                "red", "#0072B2", "#CC79A7", "#F0E442")) +
  #   ggtitle(paste0("n=", num_samp[num_iter], ", p=100")) + xlab("Index") + ylab("eigenvalues")+
  #   theme_bw()
  # gg[[num_iter]] <- print(p)
  saved_result[[num_iter]] <- tmp_mean

  cat("We are at iteration", num_iter, "\n")

}


grid.arrange(gg[[1]]$plot, gg[[2]]$plot, gg[[3]]$plot,
             gg[[4]]$plot, nrow = 4, ncol = 1)


