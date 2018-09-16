
########## FLASH on corshrink matrix  ##################

# library(flashr)
# library(softImpute)
# library(CorShrink)
# library(ggcorrplot)
#
# data("sample_by_feature_data")
# flash_out <- flash(sample_by_feature_data, Kmax = 60, tol = 0.01)
# save(flash_out, file = "output/flash_output.rda")
#

library(CorShrink)
flash_out <- get(load("../CorShrink-pages/output/flash_output.rda"))
data("sample_by_feature_data")

library(ggcorrplot)
new_data <- flash_out$EL %*% t(flash_out$EF)
cor_new_data <- cor(new_data)
rownames(cor_new_data) <- colnames(sample_by_feature_data)
colnames(cor_new_data) <- colnames(sample_by_feature_data)

image(cor_new_data)

gene_names <- as.character(read.table(file = "../data/GTEX_V6/gene_names_GTEX_V6.txt")[,1])
gene_names_1 <- as.character(sapply(gene_names, function(x) return(strsplit(x, "[.]")[[1]][1])))

person_label=read.table("../data/GTEX_V6/person_identifier_labels_with_numbers.txt");
samples_id <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,1]

samples_person <- sapply(samples_id, function(x) return(paste0(strsplit(as.character(x), "-")[[1]][1:2], collapse ="-")))
tissue_labels <- read.table(file = "../data/GTEX_V6/samples_id.txt")[,3]

unique_persons <- unique(samples_person)
unique_tissues <- unique(tissue_labels)
order_index <- c();
U <- unique_tissues
order_index <- c(order_index, grep("Brain", U))
order_index <- c(order_index, grep("fibroblast", U))
order_index <- c(order_index, grep("EBV", U))
order_index <- c(order_index, grep("Spleen", U))
order_index <- c(order_index, grep("Whole Blood", U))
order_index <- c(order_index, grep("Muscle - Skeletal", U))
order_index <- c(order_index, grep("Liver", U))
order_index <- c(order_index, grep("Pancreas", U))
order_index <- c(order_index, grep("Stomach", U))
order_index <- c(order_index, grep("Kidney - Cortex", U))
order_index <- c(order_index, grep("Adrenal Gland", U))
order_index <- c(order_index, grep("Colon - Transverse", U))
order_index <- c(order_index, grep("Small Intestine - Terminal Ileum", U))
order_index <- c(order_index, grep("Heart - Atrial Appendage", U))
order_index <- c(order_index, grep("Heart - Left Ventricle", U))
order_index <- c(order_index, grep("Minor Salivary Gland", U))
order_index <- c(order_index, grep("Skin - Sun Exposed", U))
order_index <- c(order_index, grep("Skin - Not Sun Exposed", U))
order_index <- c(order_index, grep("Lung", U))
order_index <- c(order_index, grep("Ovary", U))
order_index <- c(order_index, grep("Thyroid", U))
order_index <- c(order_index, grep("Pituitary", U))
order_index <- c(order_index, grep("Testis", U))
order_index <- c(order_index, grep("Nerve - Tibial", U))
order_index <- c(order_index, grep("Breast - Mammary Tissue", U))
order_index <- c(order_index, grep("Adipose - Visceral", U))
order_index <- c(order_index, grep("Adipose - Subcutaneous", U))
order_index <- c(order_index, grep("Artery - Coronary", U))
order_index <- c(order_index, grep("Artery - Tibial", U))
order_index <- c(order_index, grep("Artery - Aorta", U))
order_index <- c(order_index, grep("Esophagus - Mucosa", U))
order_index <- c(order_index, grep("Vagina", U))
order_index <- c(order_index, grep("Cervix - Endocervix", U))
order_index <- c(order_index, grep("Esophagus - Gastroesophageal Junction", U))
order_index <- c(order_index, grep("Colon - Sigmoid", U))
order_index <- c(order_index, grep("Esophagus - Muscularis", U))
order_index <- c(order_index, grep("Cervix - Ectocervix", U))
order_index <- c(order_index, grep("Fallopian", U))
order_index <- c(order_index, grep("Prostate", U))
order_index <- c(order_index, grep("Uterus", U))
order_index <- c(order_index, grep("Bladder", U))

ggcorrplot(cor_new_data[order_index, order_index], tl.cex = 8)
