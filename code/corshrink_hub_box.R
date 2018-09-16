

##################  Boxplots of CorShrink on hub matrices  ############################


library(CorShrink)
library(huge)
library(corpcor)
library(ggplot2)
library(gridExtra)


# hub_n_10_P_100 <- get(load("../output/hub/hub_n_10_P_100_results.rda"))
# hub_n_30_P_100 <- get(load("../output/hub/hub_n_30_P_100_results.rda"))
# hub_n_50_P_100 <- get(load("../output/hub/hub_n_50_P_100_results.rda"))
# hub_n_70_P_100 <- get(load("../output/hub/hub_n_70_P_100_results.rda"))
# hub_n_100_P_100 <- get(load("../output/hub/hub_n_100_P_100_results.rda"))
# hub_n_1000_P_100 <- get(load("../output/hub/hub_n_1000_P_100_results.rda"))

hub_n_10_P_100 <- get(load("../output/hub/hub_cmd_n_10_P_100_results.rda"))
hub_n_30_P_100 <- get(load("../output/hub/hub_cmd_n_30_P_100_results.rda"))
hub_n_50_P_100 <- get(load("../output/hub/hub_cmd_n_50_P_100_results.rda"))
hub_n_70_P_100 <- get(load("../output/hub/hub_cmd_n_70_P_100_results.rda"))
hub_n_100_P_100 <- get(load("../output/hub/hub_cmd_n_100_P_100_results.rda"))
hub_n_1000_P_100 <- get(load("../output/hub/hub_cmd_n_1000_P_100_results.rda"))

df <- data.frame("method" = c(rep("empirical", 100), rep("glasso", 100), rep("corpcor", 100), rep("corshrink", 100)),
                 "distance" = log(c(hub_n_10_P_100[,1], hub_n_10_P_100[,2],
                                    hub_n_10_P_100[,3], hub_n_10_P_100[,4])))

p <- ggplot(df, aes(method, distance)) + xlab("") + ylab("log(distance)") + ggtitle("n=10, p = 100")
p1 <- p + geom_boxplot()


df <- data.frame("method" = c(rep("empirical", 100), rep("glasso", 100), rep("corpcor", 100), rep("corshrink", 100)),
                 "distance" = log(c(hub_n_30_P_100[,1], hub_n_30_P_100[,2],
                                    hub_n_30_P_100[,3], hub_n_30_P_100[,4])))

p <- ggplot(df, aes(method, distance)) + xlab("") + ylab("log(distance)")+ ggtitle("n=30, p = 100")
p2 <- p + geom_boxplot()

df <- data.frame("method" = c(rep("empirical", 100), rep("glasso", 100), rep("corpcor", 100), rep("corshrink", 100)),
                 "distance" = log(c(hub_n_50_P_100[,1], hub_n_50_P_100[,2],
                                    hub_n_50_P_100[,3], hub_n_50_P_100[,4])))

p <- ggplot(df, aes(method, distance)) + xlab("") + ylab("log(distance)")+ ggtitle("n=50, p = 100")
p3 <- p + geom_boxplot()


df <- data.frame("method" = c(rep("empirical", 100), rep("glasso", 100), rep("corpcor", 100), rep("corshrink", 100)),
                 "distance" = log(c(hub_n_70_P_100[,1], hub_n_70_P_1000[,2],
                                    hub_n_70_P_100[,3], hub_n_70_P_100[,4])))

p <- ggplot(df, aes(method, distance)) + xlab("") + ylab("log(distance)")+ ggtitle("n=70, p = 100")
p4 <- p + geom_boxplot()

df <- data.frame("method" = c(rep("empirical", 100), rep("glasso", 100), rep("corpcor", 100), rep("corshrink", 100)),
                 "distance" = log(c(hub_n_100_P_100[,1], hub_n_100_P_100[,2],
                                    hub_n_100_P_100[,3], hub_n_100_P_100[,4])))

p <- ggplot(df, aes(method, distance)) + xlab("") + ylab("log(distance)")+ ggtitle("n=100, p = 100")
p5 <- p + geom_boxplot()

df <- data.frame("method" = c(rep("empirical", 100), rep("glasso", 100), rep("corpcor", 100), rep("corshrink", 100)),
                 "distance" = log(c(hub_n_1000_P_100[,1], hub_n_1000_P_100[,2],
                                    hub_n_1000_P_100[,3], hub_n_1000_P_100[,4])))

p <- ggplot(df, aes(method, distance)) + xlab("") + ylab("log(distance)")+ ggtitle("n=1000, p = 100")
p6 <- p + geom_boxplot()


grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3, ncol = 2)
