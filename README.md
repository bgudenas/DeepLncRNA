---
title: "README"
author: "Brian Gudenas"
date: "February 27, 2018"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, cache = TRUE)
```

## DeepLncRNA

DeepLncRNA is a deep learning model which predicts nuclear or cytosolic localization of lncRNAs. After cloning this repo and downloading the dependencies it can be utilized starting from either ensembl gene or transcript IDs or with a custom sequence.

## Dependencies
```{r , eval = FALSE}
## DeepLncRNA was built using h2o version h2o_3.18.0.11
## H2o models require you are running the same version therefore you will need to install the archived version
install.packages("https://cran.r-project.org/src/contrib/Archive/h2o/h2o_3.18.0.11.tar.gz", repos=NULL, type="source")

install.packages("dplyr")
source("https://bioconductor.org/biocLite.R")
biocLite(c("Biostrings","biomaRt"))
```


```{r , message = FALSE, warning=FALSE}
source("./src/DeepLncRNA.R")
IDs = c("ENSG00000224294", "ENSG00000278709" ) 

mart = biomaRt::useMart("ensembl","hsapiens_gene_ensembl")
trans = Gene2Transcript(IDs, mart)
seqs = GetSeq(trans$ensembl_transcript_id, mart)
df = FeatureExtract(seqs)
preds = DeepLncRNA(df, df$ensembl_transcript_id)

preds$gene = trans$ensembl_gene_id[match(rownames(preds), trans$ensembl_transcript_id)]
preds

```

