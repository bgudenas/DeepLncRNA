# create features from motif database -------------------------------------
# Brian Gudenas 7/3/2017

library(Biostrings)
library(biomaRt)

df = read.csv("./Data/DE_lncRNAs_kmertable.csv")
df$RNA = as.character(df$RNA)
# 
# mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# seq = getSequence(id = df$X, type="ensembl_transcript_id", seqType = "cdna", mart = mart)
# # seq = seq[match(df$X, seq$ensembl_transcript_id),  ]


motifs_dir = file.path("./Data/motifs/pwms_all_motifs")
mot_num = length(list.files(motifs_dir))

mat_mot = matrix(nrow = nrow(df), ncol = mot_num , data = 0)

for (i in 1:mot_num){
    fils = list.files(motifs_dir)[i]
    mot = try(read.table(paste0(motifs_dir, "/", fils), sep = "\t", row.names = 1, header = TRUE ), silent = FALSE)
    if (class(mot) != "try-error"){
        mot = t(mot)
        rownames(mot)[rownames(mot)== "U"] ="T"
        
    
        for (j in 1:nrow(df)){
        seq = str_replace_all(string = df$RNA[j], pattern = "U", replacement = "T")
        mat_mot[j,i] = countPWM(mot, seq, min.score = "80%") + countPWM(mot, reverse(seq), min.score = "80%")
        }
    }    
}

df = cbind(df, mat_mot)

source("../../../Brian/Projects/Localization/src/RNAfold.R")
out = RNAfold_out("./Data/MFE_lncRNAs")
mfe = as.numeric(as.character(out$MFE[match(df$X, out$V1)]))
df$MFE = mfe

write.csv(df, "./Data/lncRNA_featureset.csv")





