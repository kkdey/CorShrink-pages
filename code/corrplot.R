
#######  moving to corrplot  ###########

library(CorShrink)
data("sample_by_feature_data")
cor_data <- cor(sample_by_feature_data, use = "pairwise.complete.obs")

library(corrplot)
col2 <- c("blue", "white", "red")
corrplot(cor_data, diag = FALSE,
         col = colorRampPalette(col2)(200),
         tl.pos = "td", tl.cex = 0.5, tl.col = "black",
         rect.col = "white",na.label.col = "white",
         method = "color", type = "upper")
