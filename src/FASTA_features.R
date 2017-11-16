
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


Kmer_add = function(df, k){
    library(Biostrings)
    #df must contain "cdna" column
    #k is largest int of 
    kmer_mat = matrix( ncol = sum(4^seq(2, k)) , nrow = nrow(df), data =0)
    
    for (i in 1:nrow(df)){
        seq =DNAString(df$cdna[i])
        # len = df$transcript_length[i]
        index = 0
        for (j in 2:k) {
            kmers = oligonucleotideFrequency(seq, j)
            ind2 = index
            kmer_mat[i, (index + 1): (ind2 + (4^j) ) ] = as.vector(kmers)
            index = index + (4^j)
        }
    }
    
    nams =c()
    for (i in 2:k){
        nams =c(nams, names(oligonucleotideFrequency(seq, i)) )
        
    }
    colnames(kmer_mat) = nams
    df = cbind(df, kmer_mat)
    return(df)
}
dfKmer = Kmer_add(df, k = 5)

RNABP_Motif_add = function(df){
    library(Biostrings)
    
    motifs_dir = file.path("./Data/pwms_all_motifs")
    mot_num = length(list.files(motifs_dir))
    
    mat_mot = matrix(nrow = nrow(df), ncol = mot_num , data = 0)
    
    for (i in 1:mot_num){
        fils = list.files(motifs_dir)[i]
        #some motif files from CISBP are completely empty and will throw an error so i added a try-catch here
        mot = try(read.table(paste0(motifs_dir, "/", fils), sep = "\t", row.names = 1, header = TRUE ), silent = FALSE)
        if (class(mot) != "try-error"){
            mot = t(mot)
            # rownames(mot)
            # [1] "A" "C" "G" "U"
            ## Motifs are in units of RNA so we must complement to DNA (A -> T, C -> G, G -> C, U -> A)
            rownames(mot) = c("T","G","C", "A")
            mot = mot[order(rownames(mot)), ]
            #print(ncol(mot))
            
            for (j in 1:nrow(df)){
                #len = df$transcript_length[j]
                seq = df$cdna[j]
                counts = countPWM(mot, seq, min.score = "60%") 
                mat_mot[j,i] = counts
            }
        }    
    }
    
    df = cbind(df, mat_mot)
    return(df)
}
dfMotifs = RNABP_Motif_add(dfKmer)

library(h2o)
localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = 1, max_mem_size = "2G")
 # TODO: add AE feature extraction and warnings
