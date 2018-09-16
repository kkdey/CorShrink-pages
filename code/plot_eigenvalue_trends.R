

#############  eigenvalues trend plotter  ####################
library(ggplot2)
library(gridExtra)
num <- 20
eigenvalue_trends <- get(load("../output/eigenvalues_distribution/hub_eigenvalues_distribution.rda"))

num_samp <- c(10, 50, 100, 1000)

gg <- list()
for(i in 1:4){
  eigs.df <- data.frame ("x" = rep(1:num, 5),
                         "y" = eigenvalue_trends[[i]]$mean,
                         "color" = c(rep("empirical", num), rep("corpcor",  num), rep("glasso", num),
                                     rep("corshrink", num), rep("true", num)),
                         "type" = c(rep("A", 4*num), rep("B", num)))

  gg[[i]] <- ggplot(eigs.df, aes(x=x, y=y, colour=color, linetype = color)) + geom_line(lty = 1, lwd = 0.7) +
    scale_linetype_manual(values = c(rep("solid", 4), rep("dashed", 1))) +
    scale_colour_manual(values=c("#000000", "blue", "green", "gold",
                                 "red", "#0072B2", "#CC79A7", "#F0E442")) + xlab("") + ylab("")+
    theme_bw() + theme(legend.position="none")
}

grid.arrange(gg[[1]], gg[[2]], gg[[3]],
             gg[[4]], nrow = 2, ncol = 2)


eigenvalue_trends <- get(load("../output/eigenvalues_distribution/toeplitz_eigenvalues_distribution.rda"))

num_samp <- c(10, 50, 100, 1000)

for(i in 1:4){
  eigs.df <- data.frame ("x" = rep(1:num, 5),
                         "y" = eigenvalue_trends[[i]]$mean,
                         "color" = c(rep("empirical", num), rep("corpcor",  num), rep("glasso", num),
                                     rep("corshrink", num), rep("true", num)),
                         "type" = c(rep("A", 4*num), rep("B", num)))

  gg[[(4+i)]] <- ggplot(eigs.df, aes(x=x, y=y, colour=color, linetype = color)) + geom_line(lty = 1, lwd = 0.7) +
    scale_linetype_manual(values = c(rep("solid", 4), rep("dashed", 1))) +
    scale_colour_manual(values=c("#000000", "blue", "green", "gold",
                                 "red", "#0072B2", "#CC79A7", "#F0E442")) + xlab("") + ylab("")+
    theme_bw() theme(legend.position="none")
}

grid.arrange(gg[[5]], gg[[6]], gg[[7]],
             gg[[8]], nrow = 2, ncol = 2)



eigenvalue_trends <- get(load("../output/eigenvalues_distribution/banded_precision_eigenvalues_distribution.rda"))

num_samp <- c(10, 50, 100, 1000)

for(i in 1:4){
  eigs.df <- data.frame ("x" = rep(1:num, 5),
                         "y" = eigenvalue_trends[[i]]$mean,
                         "color" = c(rep("empirical", num), rep("corpcor",  num), rep("glasso", num),
                                     rep("corshrink", num), rep("true", num)),
                         "type" = c(rep("A", 4*num), rep("B", num)))

  gg[[(8+i)]] <- ggplot(eigs.df, aes(x=x, y=y, colour=color, linetype = color)) + geom_line(lty = 1, lwd = 0.7) +
    scale_linetype_manual(values = c(rep("solid", 4), rep("dashed", 1))) +
    scale_colour_manual(values=c("#000000", "blue", "green", "gold",
                                 "red", "#0072B2", "#CC79A7", "#F0E442")) + xlab("") + ylab("")+
    theme_bw() + theme(legend.position="none")
}

grid.arrange(gg[[9]], gg[[10]], gg[[11]],
             gg[[12]], nrow = 2, ncol = 2)

grid.arrange(gg[[1]], gg[[2]], gg[[3]], gg[[4]],
             gg[[5]], gg[[6]], gg[[7]], gg[[8]],
             gg[[9]], gg[[10]], gg[[11]], gg[[12]],
             nrow = 4, ncol=3, as.table=FALSE)

#############  save as 10 x 8 figure in pdf in R ##############


