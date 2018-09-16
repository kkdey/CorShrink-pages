

# if (!file.exists("cookbooks.zip")) {
#   download.file("http://archive.lib.msu.edu/dinfo/feedingamerica/cookbook_text.zip","cookbooks.zip")
# }
# unzip("cookbooks.zip",exdir="cookbooks")

## Create fake articles

# library(wordVectors)
# dir.create("../data/fake_cookbooks")
# ll <- list.files("../data/cookbooks/")
#
# num_files <- length(ll)
# for(i in 1:num_files){
#   samp <- sample(num_files, 1, replace = TRUE)
#   file.copy(paste0("../data/cookbooks/",ll[samp]), paste0("../data/fake_cookbooks/", i, ".txt"))
# }
#
# prep_word2vec(origin="../data/fake_cookbooks/",destination="../data/fake_pooled/fake_pooled.txt",lowercase=T,bundle_ngrams=1)
# model <- train_word2vec("../data/fake_pooled/fake_pooled.txt", "../data/fake_pooled/fake_pooled.bin")
#

dat <- get(load("../output/correlation_results_original_fake.rda"))

original_cors <- dat$original[lower.tri(dat$original)]
original_z <- 0.5 * log((1+original_cors)/(1-original_cors))

fake_z_mat <- matrix(0,100,length(original_z))
for(m in 1:100){
  tmp <- dat$fake[[m]][lower.tri(dat$fake[[m]])]
  fake_z_mat[m,] <- 0.5 * log((1+tmp)/(1-tmp))
}

sd_fake_z_mat <- apply(fake_z_mat, 2, sd)

library(ashr)

out <- ash(original_z, sd_fake_z_mat, mixcompdist = "normal")

ash_out <- (exp(2*out$result$PosteriorMean) - 1)/(exp(2*out$result$PosteriorMean) + 1)

corshrink_mat <- matrix(0, dim(dat$original)[1], dim(dat$original)[2])
corshrink_mat[lower.tri(corshrink_mat)] <- ash_out
corshrink_mat_2 <- corshrink_mat + t(corshrink_mat)  + diag(1, dim(dat$original)[1])

image(corshrink_mat_2)


hc <- hclust(dist(corshrink_mat_2))
hc$order
dd <- as.dendrogram(hc)
order.dendrogram(dd) ## the same :
stopifnot(hc$order == order.dendrogram(dd))

image(corshrink_mat_2[hc$order, hc$order], col=colorRampPalette(c("blue", "white", "red"))(100), zlim = c(-1,1))
image(dat$original[hc$order, hc$order], col=colorRampPalette(c("blue", "white", "red"))(100), zlim = c(-1,1))


plot(dat$original[lower.tri(dat$original)], corshrink_mat_2[lower.tri(corshrink_mat_2)])

library(ggplot2)
df <- data.frame("true_corr" = dat$original[lower.tri(dat$original)],
                 "corshrink" = corshrink_mat_2[lower.tri(corshrink_mat_2)])
p <- ggplot(df, aes(true_corr, corshrink))
p + geom_point()

idx <- which(corshrink_mat_2 < 0.3 & dat$original > 0.5, arr.ind = TRUE)
cbind(idx, apply(idx, 1, function(x) return(corshrink_mat_2[x[1], x[2]])),
      apply(idx, 1, function(x) return(dat$original[x[1], x[2]])))

rownames(corshrink_mat_2) <- rownames(dat$original)
colnames(corshrink_mat_2) <- colnames(dat$original)

cbind(which(dat$original > 0.5 &  corshrink_mat_2 < 0.22, arr.ind = TRUE)

sort(corshrink_mat_2["crisps",], decreasing = T)[1:20]
sort(dat$original["crisps",], decreasing = T)[1:20]

sort(corshrink_mat_2["lentil",], decreasing = T)[1:20]
sort(dat$original["lentil",], decreasing = T)[1:20]

sort(corshrink_mat_2["pollock",], decreasing = T)[1:20]
sort(corshrink_mat_2["haddock",], decreasing = T)[1:20]
sort(dat$original["haddock",], decreasing = T)[1:20]


