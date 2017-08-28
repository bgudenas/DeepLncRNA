# Parse RNAfold output ----------------------------------------------------

# MFE_out = "C:/Users/blg00/Google Drive/Localization/Data/MFE_lncRNAs"

RNAfold_out = function(MFE_out){
  library(stringr)
  df = read.delim(MFE_out, sep = ">", header = FALSE)
  IDs = df$V2[df$V2 != ""]
  
  
# match any number of digits following "-" with optional decimal point
  MFE = str_match(df$V1, "-\\d+\\.*\\d+")
  MFE = MFE[!is.na(MFE)]
  if (length(MFE) == length(IDs)){
    out = cbind(as.character(IDs), MFE)
  }
  return(as.data.frame(out))
}



