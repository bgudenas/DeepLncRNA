
# DeepLocal - Fasta pre-processing  ---------------------------------------
# Fasta to feature set
library(seqinr)
library(stringr)
seq = read.fasta("./Data/Test.fa", forceDNAtolower = FALSE, as.string = TRUE)
options(stringsAsFactors = FALSE)



df = matrix(nrow= length(seq), ncol =3)
colnames(df) = c( "cdna", "transcript_biotype","percentage_gene_gc_content")
for ( i in 1:length(seq)){
  df[i, 1] = as.character(seq[[i]])
  
  biotype = attr(seq[[i]], "Annot")
  df[i, 2] = str_split(biotype, "cdna:")[[1]][2]
  
  df[i, 3] = round(sum(alphabetFrequency(DNAString(df[i,1]), as.prob = TRUE)[2:3]),3)
}
  
df = as.data.frame(df)
df$percentage_gene_gc_content = as.numeric(df$percentage_gene_gc_content)
df$transcript_length = lapply(df$cdna, nchar)
  