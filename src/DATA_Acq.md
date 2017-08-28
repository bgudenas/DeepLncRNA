Data Acquisition and Processing
=============================
Author: Brian Gudenas

Date: 4/26/2017

# Description
Download Data RNA-seq data from fractionated cell lines from encodeproject.org and quantify transcript abundances
https://www.encodeproject.org/matrix/?type=Experiment

Data from ENCODE matrix portal is obtained  by using filters:
1. Homo Sapiens
2. Total RNA-seq
3. PolyA mRNA-seq
4. polyA depleted RNA-seq
5. fastq

Code
======================================
    #setup dirs
    mkdir ENCODE
    cd ENCODE
    mkdir src
    mkdir Data
    mkdir Data/Raw Data/Meta
    cd Data/Meta

	## following link was created on 7/30/17 as described above
    curl https://www.encodeproject.org/batch_download/type%3DExperiment%26assay_title%3Dtotal%2BRNA-seq%26assay_title%3DpolyA%2BmRNA%2BRNA-seq%26assay_title%3DpolyA%2Bdepleted%2BRNA-seq%26files.file_type%3Dfastq%26replicates.library.biosample.donor.organism.scientific_name%3DHomo%2Bsapiens > RAW_Links.txt
    head -n 1 RAW_Links.txt ##metadata link in first row
	# may need to paste metadata link in browser if curl doesnt work
    curl https://www.encodeproject.org/metadata/type=Experiment&replicates.library.biosample.donor.organism.scientific_name=Homo+sapiens&assay_title=total+RNA-seq&files.file_type=fastq/metadata.tsv > metadata.tsv
    cut -f 14 metadata.tsv  | sort | uniq
	##select all samples that are subcellular fractions
    awk 'BEGIN {FS="\t"} ($14 != "")' metadata.tsv  | cut -f 42 > Subcell_DL_urls
    ## Download all data
	cd ../Raw
    xargs -n 1 curl -O -L < ../Subcell_DL_urls

    ## Create file for paired-end samples mappings contains pair1 pair2 and pair1_number
    awk 'BEGIN {FS="\t"} ($14 != "")' metadata.tsv  | cut -f 1,36 > Pair_map.txt
    cd ..

Filter to remove duplicates
====================================
    R

    #Pair map contains duplicates that need to removed
    pair = read.table("./Pair_map.txt", header = 1, fill = TRUE, stringsAsFactors = FALSE)
    pairing = pair
    pair = pair[order(pair$accession), ]

    new =vector(mode = "character", length = 24)
    for ( i in 1:nrow(pair)) {
        pair = pair[pair$Paired != pair$File[i], ]
        new[i] = pair$File[i]
    	}
    new = new[!is.na(new)]

    pairs = cbind(new, pairing$Paired[match(new, pairing$File)])
    colnames(pairs) = c("Pair1","Pair2")
    pairs = pairs[rowSums(pairs=="")==0, ] ## remove any single-ended samples
    write.table(pairs, "Pair_map.txt", sep="\t", row.names =FALSE, quote=FALSE)

FastQC
=======================================
    for i in $(ls Data )
            do
            fastqc -t 8 Data/$i -o FastQC
    done


Transcript Quantification
=======================================
    ## build Salmon transcriptome
    wget ftp://ftp.ensembl.org/pub/release-88/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
    wget ftp://ftp.ensembl.org/pub/release-88/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz
    gunzip Homo_sapiens.GRCh38.ncrna.fa.gz 
    gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz

    ## Merge protein-coding and ncRNA into transcriptome file
    cat Homo_sapiens.GRCh38.cdna.all.fa Homo_sapiens.GRCh38.ncrna.fa > Hg38_Transcriptome.fa
    
    salmon -t Hg38_Transcriptome.fa -i Hg38_index
	
Download RNA motifs from CISBP-RNA and MEME-suite
========================================
    wget http://cisbp-rna.ccbr.utoronto.ca/tmp/Homo_sapiens_2017_07_07_1:58_pm.zip
	wget http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.15.tgz
    
	unzip Homo_sapiens_2017_07_07_1:58_pm.zip
    tar -xzvf motif_databases.12.15.tgz
	
	cp pwm_all_motifs ./Data/Raw
	cp /motif_databases/RNA/Ray2013_rbp_Homo_sapiens.meme ./Data/Raw
	
    
