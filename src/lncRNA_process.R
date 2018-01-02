master = read.csv("./Data/tables/Master.csv", header = TRUE, stringsAsFactors = FALSE)


annots = rtracklayer::readGFF("./Data/Annotation/gencodev27lncRNAs.gtf", version = 2L)
annots = dplyr::filter(annots, type == "transcript")
# remove trailing decimal
annots$transcript_id = lapply(strsplit(annots$transcript_id, "\\."), "[[", 1)

lnc_match = match(master$target_id, annots$transcript_id)


lncRNAs = master[!is.na(lnc_match), ]
library(dplyr)

## set up Expression matrix lncRNAs x cell_type
mat = matrix(nrow = length(unique(lncRNAs$target_id)), ncol = length(unique(lncRNAs$cell)), data = NA)
colnames(mat) = unique(lncRNAs$cell)
rownames(mat) = unique(lncRNAs$target_id)

## Loop through and add fold-changes, leave NA if not detected
for (i in 1:nrow(lncRNAs)){
  l2fc = lncRNAs$b[i]
  target = lncRNAs$target_id[i]
  celltype = lncRNAs$cell[i]
  
  mat[rownames(mat) == target, colnames(mat) == celltype] = l2fc
}
  
## truncate extremes to avoid saturation in heatmap
mat2 = mat
mat2[mat2 < -4] = -4
mat2[mat2 > 6] = 6


library(pheatmap)
pdf("./Figures/lncRNA_heatmap.pdf")
pheatmap(mat2, cluster_rows = FALSE, cluster_cols = FALSE, labels_row = "" )
dev.off()
rm(mat2)

table(rowMeans(mat, na.rm = TRUE) < 0)

wts = table(s2c$cell)/nrow(s2c)
wts = wts[order(names(wts)) ]

mat = mat[ , order(colnames(mat))]


df = data.frame(target_id = rownames(mat))
df$l2fc = NA
for (i in 1:nrow(mat)){
  
  df$l2fc[i] = weighted.mean(mat[i,], wts, na.rm = TRUE)
}



## get lncRNA cDNA sequences
library(biomaRt)
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")

seq = biomaRt::getSequence(id = df$target_id, type="ensembl_transcript_id", seqType = "cdna", mart = mart)
df$cdna = seq$cdna[match(df$target_id, seq$ensembl_transcript_id)]

df$transcript_length = as.numeric(lapply(df$cdna, nchar))
df = df[df$transcript_length >= 200, ]
saveRDS(df, "./Data/DE_lncRNAs.rds")
