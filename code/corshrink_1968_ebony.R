
########  code for ash on word2vec Ebony 1968  ##########

ll <- get(load("../output/1968_z_sdz.rda"))
library(ashr)
out <- ash(ll$z, ll$z_sd, control = list(maxiter = 500))
save(out, file = "../output/1968_ash_ebony.rda")
