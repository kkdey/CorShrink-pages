
################   Oscope data analysis   #########################

oscope_data <- read.csv("../data/GSE64016_H1andFUCCI_normalized_EC.csv", row.names = 1)
library(ccRemover)
cell_cycle_gene_indices <- gene_indexer(rownames(oscope_data), species = "human",
                                        name_type = "symbols" )

oscope_data_filtered <- oscope_data[cell_cycle_gene_indices,]
oscope_data_filtered_2 <- log(oscope_data_filtered+1)
oscope_data_filtered_2[oscope_data_filtered == 0] <- NA
#oscope_data_filtered_3 <- oscope_data_filtered_2[, -grep("H1", colnames(oscope_data_filtered_2))]
oscope_data_filtered_3 <- oscope_data_filtered_2[, grep("G2", colnames(oscope_data_filtered_2))]

presence_absence <- oscope_data_filtered_3
presence_absence[is.na(presence_absence)] <- 0
presence_absence[presence_absence > 0] <- 1
oscope_data_filtered_3 <- oscope_data_filtered_3[which(rowSums(presence_absence) > 2),]

plot(density(rowSums(presence_absence)))
nsamp_mat <- as.matrix(presence_absence) %*% t(as.matrix(presence_absence))


cor_mat <- cor(t(as.matrix(oscope_data_filtered_3)), use = "pairwise.complete.obs")
cor_mat[is.na(cor_mat)] <- 0

hc <- hclust(dist(cor_mat))
hc$order
dd <- as.dendrogram(hc)
order.dendrogram(dd) ## the same :
stopifnot(hc$order == order.dendrogram(dd))


cor2 <- cor_mat[hc$order, hc$order]
image(cor2, col=colorRampPalette(c("red", "white"))(10), zlim = c(0,1),
      xaxt = "n", yaxt = "n")

corshrink <- CorShrink::CorShrinkData(t(oscope_data_filtered_3), ash.control = list(
  mixcompdist = "halfuniform",
  control = list(maxiter=300)))



