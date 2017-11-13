
# DeepLocal - Fasta pre-processing  ---------------------------------------
# Fasta to feature set
library(seqinr)
library(stringr)
seq = read.fasta("./Data/Homo_sapiens_SPATA21_203_sequence.fa", forceDNAtolower = FALSE, as.string = TRUE)



df = matrix(nrow= length(seq), ncol =2)
colnames(df) = c( "cdna", "transcript_biotype")
for ( i in seq){
  df[i, 1] = as.character(seq[[i]])
  
  biotype = attr(seq[[i]], "Annot")
  df[i, 2] = str_split(biotype, "cdna:")[[1]][2]
}
  
  
  