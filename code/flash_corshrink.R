
########## FLASH on corshrink matrix  ##################

library(flashr)
library(softImpute)
library(CorShrink)
library(ggcorrplot)

data("sample_by_feature_data")
flash_out <- flash(sample_by_feature_data, Kmax = 100, tol = 0.001,
                   var_type = "constant")
save(flash_out, file = "output/flash_output.rda")
new_data <- flash_out$EL %*% t(flash_out$EF)
cor_new_data <- cor(new_data)

image(cor_new_data)
