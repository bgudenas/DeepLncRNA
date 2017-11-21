# Create PCA of SDA features ----------------------------------------------

load(file="./Data/SDA.RData")
dat$ID = rownames(df)
dat$biotype = df$transcript_biotype
test_compressed$ID = rownames(test_raw)
test_compressed$biotype = test_raw$transcript_biotype
train_compressed$ID = rownames(train_raw)
train_compressed$biotype = train_raw$transcript_biotype
dat = rbind(dat, test_compressed, train_compressed)
dat = dat[!duplicated(dat$ID), ]


library(ggplot2)
library(ggfortify)

dat$L2FC = NA
dat$L2FC = Super$b[match(dat$ID, rownames(Super)) ]
# - to avoid min/max extremes forcing most colors into middle we will truncate ends for viz
dat$L2FC[dat$L2FC <= -4] = -4
dat$L2FC[dat$L2FC >= 4] = 4


# Colo = rep("gray70", nrow(dat))
# Colo[dat$b >= 1] = "red4"
# Colo[dat$b <= -1] = "blue4"
PCA = prcomp(dat[ , 1:500], scale. = TRUE)
dat$alpha = FALSE
dat$alpha[!is.na(dat$L2FC)] = TRUE

pdf("./Figures/SDA_PCA.pdf")
autoplot(PCA, data = dat, colour = 'biotype', size=2, alpha = 'alpha') +
  scale_alpha_discrete( range =c(0.5, 1), guide =FALSE) +
  #scale_colour_gradient2() +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dev.off()