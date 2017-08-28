# Perform DE for subcell data ---------------------------------------------
# Brian Gudenas 6/12/2017
# import transcript quantification files produced by Salmon ---------------
library(GenomicFeatures)
library(DESeq2)
library(tximport)
library(readr)
library(stringr)


samples = read_tsv("./Data/Meta/metadata.tsv")
colnames(samples) = make.names(colnames(samples))
Pair = read.table("./Data/Meta/Pair_map.txt", sep="\t", header = TRUE)

samples = samples[!is.na(match(samples$File.accession, Pair$Pair1)), ]

samples = samples[samples$Biosample.subcellular.fraction.term.name == "nucleus" | samples$Biosample.subcellular.fraction.term.name == "cytosol", ]
## select for total RNA
table(samples$Library.made.from)
# samples = samples[samples$Library.depleted.in == "rRNA", ]

files = file.path("./Data/quants", paste0(samples$File.accession,"_quant"), "quant.sf")

txdb <- makeTxDbFromGFF("./Data/Annotation/gencode.v26.annotation.gtf.gz")
k <- keys(txdb, keytype = "GENEID")
df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene =df[ ,2:1]
# tx2gene =df[ ,c(2,2)]
# 
# head(tx2gene)
tx2gene$GENEID = unlist(lapply(str_split(tx2gene$GENEID, "\\."), "[[", 1))
head(tx2gene)

txi.salmon =tximport(files, type="salmon", tx2gene = tx2gene, reader = read_tsv)
# txi.salmon$counts = txi.salmon$counts[ rowVars(txi.salmon$counts) >= quantile(rowVars(txi.salmon$counts))[2], ]

samples$Biosample.subcellular.fraction.term.name = as.factor(samples$Biosample.subcellular.fraction.term.name)
samples$Biosample.subcellular.fraction.term.name = droplevels(samples$Biosample.subcellular.fraction.term.name)
# samples$Platform = droplevels(samples$Platform)
colnames(samples)[colnames(samples) == "Biosample.subcellular.fraction.term.name"] = "Localization"

## re-name poly-A samples which had an NA for library.depleted method
samples$Library.depleted.in[is.na(samples$Library.depleted.in)] = "rRNA"
samples$Lib_method = paste0(samples$Library.made.from, samples$Library.depleted.in)

#Check to see we have the right 3 lib types (poly-A, total RNA and poly-A-depleted)
table(samples$Lib_method)
samples$Read.length = as.factor(samples$Read.length)
## include cell type, Lib method and Platform method as confounders
dds = DESeqDataSetFromTximport(txi.salmon, samples, design = ~ Biosample.term.name + Lib_method + Platform + Localization)
dds = estimateSizeFactors(dds)
keep = rowSums(counts(dds, normalized = TRUE) >=1) >=3
dds = dds[keep, ]
dds_DE =DESeq(dds)

#load(file="./Data/DDS_trans.RData")
##contrast takes 3 args: name of variable, factor level of numerator and factor level of denominator 
res05 =results(dds_DE, alpha = 0.05, contrast = c("Localization","nucleus","cytosol")) ## 

sapply(split(colSums(counts(dds))/1e6, samples$Localization), summary)
# cytosol nucleus
# Min.      14.29   12.95
# 1st Qu.   26.52   31.53
# Median    78.55   54.00
# Mean      73.59   58.02
# 3rd Qu.  101.40   74.34
# Max.     215.50  175.80


pdf("./Figures/MA_plot.pdf")
plotMA(res05, main="MA plot: nuclear/cytoplasmic")
dev.off()

rld <- rlog(dds, blind=FALSE)
pdf("./Figures/PCA.pdf")
plotPCA(rld, intgroup=c( "Biosample.term.name"))
save.image("./Data/DEG.RData")

boxplot(assay(rld), las=2, col=c("blue","brown")[as.numeric(as.factor(samples$Localization))], names=samples$Biosample.term.name)

library("RColorBrewer")
library("pheatmap")
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(samples$Biosample.term.name, rld$Localization, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# select = order(rowMeans(counts(dds, normalized = FALSE)), decreasing = TRUE)[1:500]
# df = as.data.frame(colData(dds)[, c("Biosample.term.name","Library.depleted.in")])
# pheatmap(assay(rld)[select, ], cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = TRUE, annotation_col = df)
# 

res05= res05[!is.na(res05$padj), ] 
# res05 = res05[res05$padj <= 0.05, ]
trans = unlist(lapply(str_split(rownames(res05), "\\."), "[[", 1))
rownames(res05) = trans
         
library(biomaRt)
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map = getBM(mart = mart, attributes = c("ensembl_gene_id", "percentage_gene_gc_content","gene_biotype","chromosome_name", "external_gene_name","transcript_length", "transcript_count"),
            filters = "ensembl_gene_id", values=rownames(res05))

sig_match = match(rownames(res05),  map$ensembl_gene_id)
res05$gene_name = map$external_gene_name[sig_match]
res05$biotype = map$gene_biotype[sig_match]
res05$GC = map$percentage_gene_gc_content[sig_match]
res05$transcript_length = map$transcript_length[sig_match]
res05$transcript_count = map$transcript_count[sig_match]
lncRNA_filter = c("3prime_overlapping_ncrna", "non_coding", "antisense", "antisense_RNA", "lincRNA", "ncrna_host", "processed_transcript", "sense_intronic", "sense_overlapping","macro_lncRNA", "bidirectional_promoter_lncRNA")

lnc_test=c()
for (biotype in res05$biotype) {
    lnc_test =c(lnc_test, (sum(grepl(biotype, lncRNA_filter)) > 0 ) )
}
res05$lncRNA = lnc_test
res05 = res05[res05$lncRNA==TRUE, ]
res05 = res05[res05$transcript_length >= 200, ]

seq = getSequence(id = rownames(res05), type="ensembl_gene_id", seqType = "cdna", mart = mart)

# Convert cDNA back into RNA  --------
seq$RNA = ""
for ( i in 1:nrow(seq)){
    seq$RNA[i] = as.character(RNAString(reverseComplement(DNAString(seq$cdna[i]))))
}
seq = seq[match(rownames(res05), seq$ensembl_gene_id), ]

df = as.data.frame(res05)
df$RNA = seq$RNA

## seq: genes: 12672 to 24172



# Construct K-mer count matrix for k=1 to k =6 ----------------------------
kmer_mat = matrix( ncol = 4^6 + 4^5 + 4^4 + 4^3 + 4^2 + 4^1, nrow = nrow(df), data =0)
library(Biostrings)

for (i in 1:nrow(df)){
    s =RNAString(df$RNA[i])
    kmers = oligonucleotideFrequency(s, 6)
    kmer_mat[i,1:4096 ] = as.vector(kmers)
    
    kmers = oligonucleotideFrequency(s, 5)
    kmer_mat[i,4097:(4096 + 4^5) ] = as.vector(kmers)
    
    kmers = oligonucleotideFrequency(s, 4)
    kmer_mat[i,5121:(5120 + 4^4) ] = as.vector(kmers)
    
    kmers = oligonucleotideFrequency(s, 3)
    kmer_mat[i,5377:(5376 + 4^3) ] = as.vector(kmers)
    
    kmers = oligonucleotideFrequency(s, 2)
    kmer_mat[i,5441:(5440 + 4^2) ] = as.vector(kmers)
    
    kmers = oligonucleotideFrequency(s, 1)
    kmer_mat[i,5457:(5456 + 4^1) ] = as.vector(kmers)
}
colnames(kmer_mat) = c(names(oligonucleotideFrequency(s, 6)), names(oligonucleotideFrequency(s, 5)), names(oligonucleotideFrequency(s, 4)), names(oligonucleotideFrequency(s, 3)), names(oligonucleotideFrequency(s, 2)), names(oligonucleotideFrequency(s, 1)) )


df = cbind(df, as.data.frame(kmer_mat))
df = df[ , colnames(df) != "lncRNA" & colnames(df) != "ensembl_gene_id" ]

df = df[!is.na(df$log2FoldChange), ]


# need to truncate sequence length so for the RNAfold software  (  < 32,0000) ------------
for (i in 1:nrow(df)){
  write(paste0(">", rownames(df)[i], "\n", as.character(str_sub(df$RNA[i], start = 1, end = 32000), "\n" ), "./Data/Trans.fa", append = TRUE, sep = "")
}

write.csv( df , "./Data/DE_lncRNAs_kmertable.csv", row.names = TRUE)
