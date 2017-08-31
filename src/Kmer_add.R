library(Biostrings)

Kmer_add = function(df, k){
  #df must contain "cdna" column
  #k is largest int of 
  kmer_mat = matrix( ncol = sum(4^seq(1, k)) , nrow = nrow(df), data =0)
  
  for (i in 1:nrow(df)){
    seq =DNAString(df$cdna[i])
    index = 0
    for (j in 1:k) {
    kmers = oligonucleotideFrequency(seq, j)
    ind2 = index
    if (j == 1){ ind2 = 0}
    
    kmer_mat[i, (index + 1): (ind2 + (4^j) ) ] = as.vector(kmers)
    index = index + (4^j)
    }
  }
  
  nams =c()
  for (i in 1:k){
    nams =c(nams, names(oligonucleotideFrequency(seq, i)) )
  
  }
  colnames(kmer_mat) = nams
  df = cbind(df, kmer_mat)
  return(df)
}

# lncRNAs2 = Kmer_add(lncRNAs, 6)
    