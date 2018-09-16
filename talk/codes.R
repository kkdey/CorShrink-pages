
ash_words <- c("ash", "flash", "mash", "smash", "vash", "fash", "succotash",
               "truncash", "cash", "dash", "mn-mash")

library(CorShrink)
data("sample_by_feature_data")

donor_by_tissue_PLIN1 <- sample_by_feature_data

donor_by_tissue_PLIN1[1:5,1:5]

length(which(is.na(donor_by_tissue_PLIN1)))/length(donor_by_tissue_PLIN1)

library(ggcorrplot)

order_indices <- get(load("../output/order_index.rda"))

ggcorrplot(cor(donor_by_tissue_PLIN1, use = "pairwise.complete.obs")[order_indices, order_indices]) + coord_fixed(ratio=1)

out <- CorShrinkData(donor_by_tissue_PLIN1)
ggcorrplot(out$ash_cor_PD[order_indices, order_indices]) + coord_fixed(ratio=1)


library(Matrix)
nearPD(rbind(c(1,0,0,0), c(0,1,0,0), c(0,0,0,0), c(0,0,0,0)))
