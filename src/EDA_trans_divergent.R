
# Explore if isoforms diverge in Localization -----------------------------
load(file="./Data/DEG_trans.RData")

txdb <- makeTxDbFromGFF("./Data/Annotation/gencode.v26.annotation.gtf.gz")
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene =df[ ,c(1,2)]

df = df[ , 1:6]
tx2gene$TXNAME = unlist(lapply(str_split(tx2gene$TXNAME, "\\."), "[[", 1))


genes = tx2gene$GENEID[match(df$ensembl_transcript_id, tx2gene$TXNAME)]
df$genes = genes

div = 0
for ( i in unique(df$genes)){
    trans = df[df$genes == i, ]
    if (nrow(trans) > 1){
        if (min(trans$log2FoldChange) < 0 & max(trans$log2FoldChange) > 1  ){
            div = div + 1
        }
    }
}
        

# 938 / 12449 (7.53%) LncRNA genes show transcripts with divergent localizations