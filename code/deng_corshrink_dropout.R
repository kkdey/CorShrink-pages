

########      Deng et al CorShrink on dropouts   #####################

library(devtools)
#install_github('kkdey/singleCellRNASeqMouseDeng2014')
library(singleCellRNASeqMouseDeng2014)
deng.counts <- Biobase::exprs(Deng2014MouseESC)
deng.meta_data <- Biobase::pData(Deng2014MouseESC)
deng.gene_names <- rownames(deng.counts)

feature_depth <- colSums(deng.counts)
boxplot(feature_depth ~ as.factor(deng.meta_data[,1]))

presence_absence <- deng.counts
presence_absence[presence_absence > 0] = 1
boxplot(colSums(presence_absence) ~ as.factor(deng.meta_data[,1]))

plot(density(rowSums(presence_absence)))
plot(density(rowSums(deng.counts)))

deng.counts_filtered <- deng.counts[which(rowSums(deng.counts) >  2),]
presence_absence_filtered <- presence_absence[which(rowSums(presence_absence) >  2),]

plot(density(rowSums(deng.counts_filtered)))
plot(density(rowSums(presence_absence_filtered)))

TF_genes <- as.character(read.table("../data/all_genes_deng_tf.txt")[,1])

deng.counts_filtered_2 <- deng.counts_filtered[match( TF_genes, rownames(deng.counts_filtered)),]
library(limma)
voom_deng.counts_filtered_2 <- voom(deng.counts_filtered_2)$E
voom_deng.counts_filtered_2[deng.counts_filtered_2==0] <- NA

cor_mat <- cor(t(voom_deng.counts_filtered_2), use = "pairwise.complete.obs")

write.csv(deng.counts_filtered_2, file = "../data/deng_filtered.csv")


voom_deng.counts_filtered_3 <- DrImpute(log(deng.counts_filtered_2+1))

cor_mat_2 <- cor(t(voom_deng.counts_filtered_3), use = "pairwise.complete.obs")

voom_deng.counts_filtered_4 <- voom_deng.counts_filtered_3
voom_deng.counts_filtered_4[voom_deng.counts_filtered_4==0] <- NA

out <- CorShrink::CorShrinkData(t(voom_deng.counts_filtered_4), sd_boot = FALSE,
                                ash.control = list(control = list(maxiter=100)))

save(out, file = "../output/corshrink_deng_TF_DrImpute.rda")





deng.counts_filtered_3 <- scimpute(count_path = "../data/deng_filtered.csv",
                                   infile = "csv",           # format of input file
                                   outfile = "csv",          # format of output file
                                   out_dir = "../data/",
                                   drop_thre = 0.5,          # threshold set on dropout probability
                                   Kcluster = 2)




library(ggcorrplot)

#out <- CorShrink::CorShrinkData(t(voom_deng.counts_filtered_2), sd_boot = FALSE,
#                                ash.control = list(control = list(maxiter=1000)))

#save(out, file = "../output/corshrink_deng_TF.rda")

out <- get(load(file = "../output/corshrink_deng_TF.rda"))
presence_absence <- t(voom_deng.counts_filtered_2)
presence_absence[presence_absence > 0] = 1
presence_absence[is.na(presence_absence)] = 0
nsamp_mat <- t(presence_absence) %*% presence_absence

df <- data.frame("empirical" = as.vector(cor_mat[lower.tri(cor_mat)]),
                 "corshrink" = as.vector(out$ash_cor_PD[lower.tri(out$ash_cor_PD)]),
                 "color" = as.vector(nsamp_mat[lower.tri(nsamp_mat)]))

p <- ggplot(df, aes(empirical, corshrink)) +
  geom_point(aes(colour = color)) +
  scale_x_continuous(limits = c(-1, 1)) +
  scale_y_continuous(limits = c(-1, 1))

#ggsave(filename="../figures/deng_TF.pdf", plot=p)


hc <- hclust(dist(out$ash_cor_PD))
hc$order
dd <- as.dendrogram(hc)
order.dendrogram(dd) ## the same :
stopifnot(hc$order == order.dendrogram(dd))

cor2 <- cor_mat[hc$order, hc$order]
image(cor2, col=colorRampPalette(c("blue", "white", "red"))(50), zlim = c(-1,1))

cor3 <- out$ash_cor_PD[hc$order, hc$order]

col=c(rep(rgb(0,1,0), 1),rev(rgb(seq(1,0,length=1000),1,seq(1,0,length=1000))),
      rgb(1,seq(1,0,length=1000),seq(1,0,length=1000)), rep(rgb(1,0,0),1))

image(cor3, col=colorRampPalette(c("blue", "white", "red"))(100), zlim = c(-0.7,0.7),
      xaxt = "n", yaxt = "n")







t3  <- system.time(image(cor_mat))

hmcols<-colorRampPalette(c("red","white","blue"))(256)
heatmap(cor_mat, col = hmcols)

library(gplots)
tmp <- heatmap(cor_mat, col = colorRampPalette(c("blue", "white", "red"))(20), Rowv = NA, Colv = NA,
                 trace = "none", cexRow = 0.3, cexCol = 0.3, symm=TRUE)
tmp$rowInd



