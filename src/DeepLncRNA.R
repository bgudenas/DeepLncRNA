# DeepLncRNA.R -------------------------------------------------------------

Gene2Transcript = function(Ensembl_gene_IDs , mart){
    # Extract all transcripts for a list of ensembl gene IDs
  #mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    map = biomaRt::getBM(mart = mart, attributes = c("ensembl_gene_id", "ensembl_transcript_id","transcript_biotype"), filters = "ensembl_gene_id", values = Ensembl_gene_IDs)
    return(map)
}


GetSeq = function(Ensembl_transcript_IDs, mart, versioned=FALSE){
    ## Uses transcript_Ids to extract, biotype, chromosome, cDNA sequence, transcript_length, GC %
    library(Biostrings)
    #mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    if(versioned == TRUE) {
      map = biomaRt::getBM(mart = mart, attributes = c( "ensembl_transcript_id_version","transcript_biotype", "chromosome_name"), filters = "ensembl_transcript_id_version", values = Ensembl_transcript_IDs)
      colnames(map)[1] <- "ensembl_transcript_id"
      seq = biomaRt::getSequence(id = map$ensembl_transcript_id, type="ensembl_transcript_id_version", seqType = "cdna", mart = mart)
    } else {
      map = biomaRt::getBM(mart = mart, attributes = c( "ensembl_transcript_id","transcript_biotype", "chromosome_name"), filters = "ensembl_transcript_id", values = Ensembl_transcript_IDs)
      seq = biomaRt::getSequence(id = map$ensembl_transcript_id, type="ensembl_transcript_id", seqType = "cdna", mart = mart)
    }
    map$cdna = seq$cdna[match(map$ensembl_transcript_id, seq$ensembl_transcript_id)]
    
    map$transcript_length = as.numeric(lapply(map$cdna, nchar))
    
    for (i in 1:nrow(map)){
        map$GC_content[i] =  sum(oligonucleotideFrequency(DNAString(map$cdna[i]), 1)[c(2,3)])/ sum(oligonucleotideFrequency(DNAString(map$cdna[i]), 1))
    }
    return(map)
}


Kmer_add = function(df, k){
    library(Biostrings)
    #df must contain "cdna" column
    #k is largest int of k-mer
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


RNABP_Motif_add = function(df){
    library(Biostrings)
    
    motifs_dir = file.path("./Data/Motifs/pwms_all_motifs")
    mot_num = length(list.files(motifs_dir))
    
    mat_mot = matrix(nrow = nrow(df), ncol = mot_num , data = 0)
    
    for (i in 1:mot_num){
        fils = list.files(motifs_dir)[i]
        #some motif files from CISBP are empty (RFAM) and will throw an error 
        mot = try(read.table(paste0(motifs_dir, "/", fils), sep = "\t", row.names = 1, header = TRUE ), silent = TRUE)
        if (class(mot) != "try-error"){
            mot = t(mot)
            rownames(mot)[4] = "T"
            
            for (j in 1:nrow(df)){
                #len = df$transcript_length[j]
                seq = df$cdna[j]
                counts = countPWM(mot, seq, min.score = "80%") 
                mat_mot[j,i] = counts
            }
        }    
    }
    
    df = cbind(df, mat_mot)
    return(df)
}


FeatureExtract = function(map){
    library(dplyr)
    ## Use cDNA sequence to compute kmer frequencies and RNA binding motifs
    df = Kmer_add(map, k = 5)
    df = RNABP_Motif_add(df)
    #########################
    # Add biotype and chromosome feature
    df$lincRNA = 0
    df$lincRNA[df$transcript_biotype == "lincRNA"] = 1
    
    df$antisense = 0
    df$antisense[df$transcript_biotype == "antisense"] = 1
    
    df$sense = 0
    df$sense[grepl("^sense_", df$transcript_biotype)] = 1
    
    df = df %>% dplyr::select(-transcript_biotype)
    df$chromosome_name = as.factor(df$chromosome_name)
    
    chroms = data.frame(matrix(nrow=nrow(df), ncol = 24, data = 0))
    colnames(chroms) = paste0("chromosome_name", c(1:22, "X","Y"))
    for ( i in 1:nrow(df)){
        chroms[i, colnames(chroms) == paste0("chromosome_name", df$chromosome_name[i]) ] = 1
    }

    
    df = cbind(df, chroms)
    df = df %>% dplyr::select( - chromosome_name)
    
    return(df)
}


DeepLncRNA = function(df, IDs){
  library(h2o)
  localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = 1, max_mem_size = "1G")
  DNN = h2o.loadModel("./Data/Models/DeepLncRNA")
  pred_hex = as.h2o(df)
  vals = as.data.frame(predict(DNN, pred_hex, "probs"))
  rownames(vals) = IDs
  h2o.shutdown(FALSE)
  return(vals)
}

##Usage example
# #Usage
# trans = Gene2Transcript(c("ENSG00000229807", "ENSG00000228630"))
# lncRNAs = GetSeq(trans$ensembl_transcript_id)
# df = FeatureExtract(lncRNAs)
# preds = DeepLncRNA(df)    

