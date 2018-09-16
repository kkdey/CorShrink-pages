

############   CorShrink on Chu et al 2016  data  #######################

scData <- read.csv("../data/Chu2016/GSE75748_sc_cell_type_ec.csv", row.names = 1)
sc_timecourse <- read.csv("../data/Chu2016/GSE75748_sc_time_course_ec.csv", row.names = 1)

bulkData <- read.csv("../data/Chu2016/GSE75748_bulk_cell_type_ec.csv", row.names = 1)
bulk_timecourse <- read.csv("../data/Chu2016/GSE75748_bulk_time_course_ec.csv", row.names = 1)

cor_cells_timecourse <- cor(sc_timecourse)
library(ggcorrplot)
ggcorrplot(cor_cells_timecourse)

cor_bulk_timecourse <- cor(bulk_timecourse)
ggcorrplot(cor_bulk_timecourse)

sc_timecourse_filtered <- sc_timecourse[,93:dim(sc_timecourse)[2]]
cor_cells_timecourse <- cor(sc_timecourse_filtered)
image(cor_cells_timecourse)

time_span <- sapply(colnames(sc_timecourse_filtered), function(x) return(strsplit(x, "_")[[1]][1]))
boxplot(cor_cells_timecourse ~ as.factor(time_span))

num_zeros <- apply(sc_timecourse_filtered, 2, function(x) return (length(which(x==0))))
boxplot(num_zeros ~ as.factor(time_span))
