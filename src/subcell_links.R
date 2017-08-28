DL_links = read.table("./RNA-seq_URL_list.txt")
subcell_files = read.table("./Meta_subcell_filt.tsv")

keep_vec =c()

for (i in subcell_files$File.accession){
    keep_vec =c(keep_vec, as.character(DL_links[grepl(x = DL_links[,1], pattern = i), ] ))
}
    
keep_vec
write(keep_vec, file="subcell_DL_links.txt")
