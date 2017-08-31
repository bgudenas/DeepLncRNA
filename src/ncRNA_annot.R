source("https://bioconductor.org/biocLite.R")
biocLite("seqinr")
library(seqinr)
library(stringr)
library(biomaRt)

ncRNA = read.fasta("./Data/Homo_sapiens.GRCh38.ncrna.fa")
Annot = getAnnot.list(ncRNA)

out =c()

for (i in Annot){
  ID = strsplit(i[1][[1]], split = "\\.")[[1]][1]
  ID = str_sub(ID, 2, nchar(ID))
  biotype = strsplit(i[1][[1]], split = "transcript_biotype:")[[1]][2]
  biotype = strsplit(biotype, "gene_symbol:")[[1]][1]
  biotype = str_trim(biotype)
  
  lncRNA_filter = c("3prime_overlapping_ncrna", "antisense","antisense_RNA", "lincRNA","ncrna_host","processed_transcript", "sense_intronic" , "sense_overlapping", "bidirectional_lncRNA", "non_coding", "macro_lncRNA")
  
  if (biotype %in% lncRNA_filter){
    out =c(out, ID)
  }
}
    
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")


seq = getSequence(id = map$ensembl_transcript_id, type="ensembl_transcript_id", seqType = "cdna", mart = mart)

  


  