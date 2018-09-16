

################   Box plot on banded sparse matrices #########################

##################  Boxplots of CorShrink on banded_precision matrices  ############################


library(CorShrink)
library(huge)
library(corpcor)
library(ggplot2)
library(gridExtra)


banded_precision_n_10_P_100 <- get(load("../output/banded_precision/banded_precision_n_10_P_100_results.rda"))
banded_precision_n_30_P_100 <- get(load("../output/banded_precision/banded_precision_n_30_P_100_results.rda"))
banded_precision_n_50_P_100 <- get(load("../output/banded_precision/banded_precision_n_50_P_100_results.rda"))
banded_precision_n_70_P_100 <- get(load("../output/banded_precision/banded_precision_n_70_P_100_results.rda"))
banded_precision_n_100_P_100 <- get(load("../output/banded_precision/banded_precision_n_100_P_100_results.rda"))
banded_precision_n_1000_P_100 <- get(load("../output/banded_precision/banded_precision_n_1000_P_100_results.rda"))



df <- data.frame("method" = c(rep("empirical", 100), rep("glasso", 100), rep("corpcor", 100), rep("corshrink", 100)),
                 "distance" = log(c(banded_precision_n_10_P_100[,1], banded_precision_n_10_P_100[,2],
                                    banded_precision_n_10_P_100[,3], banded_precision_n_10_P_100[,4])))

p <- ggplot(df, aes(method, distance)) + xlab("") + ylab("log(distance)") + ggtitle("n=10, p = 100")
p1 <- p + geom_boxplot()


df <- data.frame("method" = c(rep("empirical", 100), rep("glasso", 100), rep("corpcor", 100), rep("corshrink", 100)),
                 "distance" = log(c(banded_precision_n_30_P_100[,1], banded_precision_n_30_P_100[,2],
                                    banded_precision_n_30_P_100[,3], banded_precision_n_30_P_100[,4])))

p <- ggplot(df, aes(method, distance)) + xlab("") + ylab("log(distance)")+ ggtitle("n=30, p = 100")
p2 <- p + geom_boxplot()

df <- data.frame("method" = c(rep("empirical", 100), rep("glasso", 100), rep("corpcor", 100), rep("corshrink", 100)),
                 "distance" = log(c(banded_precision_n_50_P_100[,1], banded_precision_n_50_P_100[,2],
                                    banded_precision_n_50_P_100[,3], banded_precision_n_50_P_100[,4])))

p <- ggplot(df, aes(method, distance)) + xlab("") + ylab("log(distance)")+ ggtitle("n=50, p = 100")
p3 <- p + geom_boxplot()


df <- data.frame("method" = c(rep("empirical", 100), rep("glasso", 100), rep("corpcor", 100), rep("corshrink", 100)),
                 "distance" = log(c(banded_precision_n_70_P_100[,1], banded_precision_n_70_P_100[,2],
                                    banded_precision_n_70_P_100[,3], banded_precision_n_70_P_100[,4])))

p <- ggplot(df, aes(method, distance)) + xlab("") + ylab("log(distance)")+ ggtitle("n=70, p = 100")
p4 <- p + geom_boxplot()

df <- data.frame("method" = c(rep("empirical", 100), rep("glasso", 100), rep("corpcor", 100), rep("corshrink", 100)),
                 "distance" = log(c(banded_precision_n_100_P_100[,1], banded_precision_n_100_P_100[,2],
                                    banded_precision_n_100_P_100[,3], banded_precision_n_100_P_100[,4])))

p <- ggplot(df, aes(method, distance)) + xlab("") + ylab("log(distance)")+ ggtitle("n=100, p = 100")
p5 <- p + geom_boxplot()

df <- data.frame("method" = c(rep("empirical", 100), rep("glasso", 100), rep("corpcor", 100), rep("corshrink", 100)),
                 "distance" = log(c(banded_precision_n_1000_P_100[,1], banded_precision_n_1000_P_100[,2],
                                    banded_precision_n_1000_P_100[,3], banded_precision_n_1000_P_100[,4])))

p <- ggplot(df, aes(method, distance)) + xlab("") + ylab("log(distance)")+ ggtitle("n=1000, p = 100")
p6 <- p + geom_boxplot()


grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3, ncol = 2)
