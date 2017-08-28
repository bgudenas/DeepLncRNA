# Deep Learning for lncRNA localization == AUTOENCODER
# 2/27/16

library(h2o)
library(genefilter)
set.seed(54321)

# Import data -------------------------------------------------------------
df = read.csv("./Data/lncRNA_featureset.csv", row.names = 1)
df = df[!is.na(df$log2FoldChange), ]
filt = -c(1,2,4,5,6,7,8,13)
df = df[ , filt]


Loc = cbind(local = rep("", nrow(df)), df)
Loc$local = as.character(Loc$local)
Loc$local[Loc$log2FoldChange <= 0] = "Cytosol"
Loc$local[Loc$log2FoldChange >= 2] = "Nuclear"
Loc = Loc[ ,-2]

col_vars = rowVars(t(Loc[ , 3:ncol(Loc)]))
Loc_cols = Loc[ , 3:ncol(Loc)]
Loc_cols = Loc_cols[, col_vars > quantile(col_vars, .2)]

Loc = cbind(Loc[ , 1:2], Loc_cols)

# Checl=k factor var for single levels ------------------------------------
table(Loc$biotype)
Loc$biotype[Loc$biotype == "macro_lncRNA"] = "lincRNA"
Loc$biotype[Loc$biotype == "non_coding"] = "lincRNA"
Loc$biotype = droplevels(Loc$biotype)
Loc$local = as.factor(Loc$local)


Loc = as.h2o(Loc)

# Loc = Loc[Loc$local != "", ]
Loc_un = Loc[Loc$local == "", ]

localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = -1, max_mem_size = "14G")


AE = h2o.deeplearning( x = 2:ncol(Loc_un), #features
                       training_frame = Loc_un,
                       activation = "Tanh",
                       input_dropout_ratio = 0.1,
                       seed = 54321,
                       hidden = c(500, 500),
                       ignore_const_cols = FALSE,
                       epochs = 100,
                       autoencoder = TRUE,
                       model_id="ae_model")

h2o.saveModel(object = AE, path = "./Data/", force = TRUE)
h2o

Loc =  Loc[ Loc$local != "", ]
splits = h2o.splitFrame(Loc, ratios =  c(0.7), seed = 54321)
train = splits[[1]]
test = splits[[2]]

H2O = h2o.deeplearning( x = 2:ncol(train), #features
                       y=1,
                       training_frame = train,
                       validation_frame = test,
                       input_dropout_ratio = 0.2,
                       l2 = 1e-5,
                       seed = 54321,
                       hidden = c(500),
                       epochs = 500,
                       ignore_const_cols = FALSE,
                       balance_classes = TRUE,
                       pretrained_autoencoder = "ae_model",
                       model_id = "Sup_AE"
                      )


