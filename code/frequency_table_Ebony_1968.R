

##################  Word counts and word presence-absence magazines  #########################


files <- list.files("../data/Ebony_1968/")
Sys.setlocale(locale="C")

all_names <- c()
for(m in 1:length(files)){
  sentences<-scan(paste0("../data/Ebony_1968/", files[m]),"character",sep="\n");
  #Replace full stop and comma
  sentences<-gsub("\\.","",sentences)
  sentences<-gsub("\\,","",sentences)
  #Split sentence
  words<-strsplit(sentences," ")
  #Calculate word frequencies
  words.freq<-table(unlist(words));
  all_names <- union(all_names , names(words.freq))
}

tab <- matrix(0, length(files), length(all_names))
colnames(tab) <- all_names
rownames(tab) <- files

for(m in 1:length(files)){
  sentences<-scan(paste0("../data/Ebony_1968/", files[m]),"character",sep="\n");
  #Replace full stop and comma
  sentences<-gsub("\\,","",sentences)
  #Split sentence
  words<-strsplit(sentences," ")
  #Calculate word frequencies
  words.freq<-table(unlist(words));
  tab[m, match(names(words.freq), colnames(tab))] <- as.numeric(words.freq)
}

save(tab, file = "../output/word2vec_Ebony/word_frequencies.rda")

