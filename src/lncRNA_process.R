master = read.csv("./Data/tables/Master.csv", header = TRUE)



lncRNAs = rtracklayer::readGFF("./Data/Annotation/gencodev27lncRNAs.gtf", version = 2L)
lncRNAs = dplyr::filter(lncRNAs, type == "transcript")
# remove trailing decimal
lncRNAs$transcript_id = lapply(strsplit(lncRNAs$transcript_id, "\\."), "[[", 1)

lnc_match = match(master$target_id, lncRNAs$transcript_id)


rawDElncRNAs = master[!is.na(lnc_match), ]
library(dplyr)

DElncRNAs =  dplyr::select(rawDElncRNAs, target_id, b, cell)  %>% tidyr::spread( key  = target_id, b)
