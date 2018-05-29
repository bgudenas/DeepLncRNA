
# Parse LncRNAs from LncLocator -------------------------------------------
# source("https://bioconductor.org/biocLite.R")
# biocLite("seqinr")

library(stringr)
library(seqinr)
library(Biostrings)
cyto = read.fasta("./Data/LncLocator/Cytosol", as.string = TRUE)
nuc = read.fasta("./Data/LncLocator/Nucleus", as.string = TRUE)
seqs =c( nuc, cyto)
nucnames = unlist(getAnnot(nuc))
cytonames = unlist(getAnnot(cyto))

allnames = c(nucnames, cytonames)
NRs = grep( "NR_*", allnames)
refseq = str_extract( allnames[NRs], "NR_\\d*") ##get Refseq IDs
refseq2 = str_extract( allnames, "NM_\\d*") ##get Refseq IDs
refseq2 = refseq2[!is.na(refseq2)]

# mice = str_extract( allnames, "XR_\\d*") 
# mice = mice[!is.na(mice)] 
# genmice = str_extract( allnames, "AK\\d*") 
# genmice = genmice[!is.na(genmice)]
ensem= str_extract( allnames, "ENST\\d*") ##get Ensembl IDs
ensem = ensem[!is.na(ensem)]

fullnames = c(refseq, refseq2, ensem)

## Get transcript biotype and chromosome if available
mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map = biomaRt::getBM(mart = mart, attributes = c( "gene_biotype","refseq_ncrna", "chromosome_name"), filters = "refseq_ncrna", values = refseq)
mapNM = biomaRt::getBM(mart = mart, attributes = c( "gene_biotype","refseq_ncrna", "chromosome_name"), filters = "refseq_mrna", values = refseq2)
map2 = biomaRt::getBM(mart = mart, attributes = c( "gene_biotype","ensembl_transcript_id" ,"chromosome_name"), filters = "ensembl_transcript_id", values = ensem)
colnames(map)= c("gene_biotype","ID", "chromosome_name")
colnames(mapNM)= c("gene_biotype","ID", "chromosome_name")
colnames(map2)= c("gene_biotype","ID", "chromosome_name")
maps = rbind(map, mapNM, map2)
maps = maps[nchar(maps$chromosome_name) <3, ]
maps = maps[!is.na(maps$ID), ]

#Loop through to make DF
fullnames = data.frame("IDs" = c(refseq, refseq2, ensem), transcript_biotype = "", chromosome_name = "",GC_content= 0, cdna = "", stringsAsFactors = FALSE)
for (i in 1:nrow(fullnames)){
  
    bio = maps$gene_biotype[maps$ID == fullnames$IDs[i]]
        if (length(bio) != 0 ){ 
            fullnames$transcript_biotype[i] = bio 
        } else {fullnames$transcript_biotype[i] = NA}
    
    chrom = maps$chromosome_name[maps$ID == fullnames$IDs[i]]
        if (length(chrom) != 0 ){
            fullnames$chromosome[i] =  chrom 
            } else {fullnames$transcript_biotype[i] = NA}

    fullnames$cdna[i] = toupper(as.character(seqs[grepl(fullnames$IDs[i], allnames)]))
    fullnames$GC_content[i] = sum(oligonucleotideFrequency(DNAString(fullnames$cdna[i]), 1)[c(2,3)])/ sum(oligonucleotideFrequency(DNAString(fullnames$cdna[i]), 1))
}
fullnames$transcript_length = as.numeric(lapply(fullnames$cdna, nchar))
fullnames = fullnames[ ,colnames(fullnames)!= "chromosome"]


allseqs = data.frame("IDs" = c(allnames), transcript_biotype = "", chromosome_name = "",GC_content= 0, cdna = "", stringsAsFactors = FALSE)
for (i in 1:nrow(allseqs)) {
    allseqs$IDs[i] = allnames[i]
    allseqs$cdna[i] = toupper(as.character(seqs[allseqs$IDs[i] == allnames]))
    allseqs$GC_content[i] = sum(oligonucleotideFrequency(DNAString(allseqs$cdna[i]), 1)[c(2,3)])/ sum(oligonucleotideFrequency(DNAString(allseqs$cdna[i]), 1))
}
allseqs$transcript_length = as.numeric(lapply(allseqs$cdna, nchar))

## FIlter out seqs from allseqs for which we got genomic features
for(i in fullnames$IDs){
    allseqs = allseqs[!grepl(i, allseqs$IDs), ]
}

master = rbind(fullnames, allseqs)

master$Loc = ""
for (i in 1:nrow(master)){
    ID = master$IDs[i]
    if (sum(grepl(ID, cytonames)) == 1 ){master$Loc[i] = "Cytosol"
        }else if (sum(grepl(ID, nucnames)) == 1 ){master$Loc[i] = "Nuclear"
            }else if (sum(ID == nucnames)==1) { master$Loc[i] = "Nuclear"
                }else if (sum(ID == cytonames)==1) { master$Loc[i] = "Cytosol" }
}
    

source("./src/DeepLncRNA.R")
lncs = FeatureExtract(master)

## Remove any lncRNAs present in training or validation data
filt = is.na(match(master$IDs, c(rownames(res$train), rownames(res$validate) )) )
lncs = lncs[filt, ]
preds = DeepLncRNA(lncs, lncs$IDs)

caret::confusionMatrix(data = preds$predict, reference = as.factor(master$Loc[filt]), positive = "Nuclear")
