


load(file="./Data/SDA.RData")


d = density(Super$b)
plot(d, main = "Distribution of lncRNA log2 fold changes")
polygon(d, col = "blue")
dat.pca = prcomp(dat, scale. = TRUE, center = TRUE)

library(ggplot2)
library(ggfortify)

dat$L2FC = NA
dat$L2FC = Super$b[match(rownames(df), rownames(Super)) ]
# - to avoid min/max extremes forcing most colors into middle we will truncate ends for viz
dat$L2FC[dat$L2FC <= -3] = -3
dat$L2FC[dat$L2FC >= 3] = 3


# Colo = rep("gray70", nrow(dat))
# Colo[dat$b >= 1] = "red4"
# Colo[dat$b <= -1] = "blue4"
df = prcomp(dat[ , -ncol(dat)], scale. = TRUE)
dat$alpha = FALSE
dat$alpha[!is.na(dat$L2FC)] = TRUE

pdf("./Figures/SDA_PCA.pdf")
autoplot(df, data = dat, colour = 'L2FC', size=2, alpha = 'alpha') +
  scale_alpha_discrete( range =c(0.3, 1), guide =FALSE) +
  scale_colour_gradient2() +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dev.off()