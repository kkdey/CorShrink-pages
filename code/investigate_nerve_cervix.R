
######### investigating Nerve Tibial and Cervix Endovervix  ###############

corr_tissues <- get(load("../output/cor_tissues_non_ash_voom_pearson.rda"))
corr_tissues[25,39,]

out <- ashr::ash(corr_tissues[52, 39, ], rep(1, dim(corr_tissues)[3]),
                 mixcompdist = "halfuniform", control = list(maxiter = 1000))
order_index <- get(load("../output/order_index.rda"))
