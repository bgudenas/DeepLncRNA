# Deep Learning for lncRNA localization TRANSCRIPTS
# 2/27/16


library(h2o)
library(genefilter)
set.seed(54321)

# Import data -------------------------------------------------------------
df = read.csv("./Data/DE_lncRNAs_kmertable_Trans.csv", row.names = 1)
filt = -c(6,7)
df = df[ , filt]

feature_filt = rowVars(t(df[ ,6:ncol(df)]))
kmers = df[ ,6:ncol(df)]
kmers = kmers[ , feature_filt >= quantile(feature_filt)[2]]
dat = cbind(df[, 1:5], kmers)
dat$log2FoldChange = round(dat$log2FoldChange, 4)
# Localization = rep("Equal", nrow(df))
# Localization[df$log2FoldChange < 0 ] = "Cytosol"
# Localization[df$log2FoldChange > 2 ] = "Nuclear"
# df$log2FoldChange = as.factor(Localization)
rm(kmers,df)


# Start H2O cluster -------------------------------------------------------
localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = -1, max_mem_size = "14G")

dat_h2o = as.h2o(dat)
splits = h2o.splitFrame(dat_h2o, ratios =  c(0.7, 0.15), seed = 54321)
train = splits[[1]]
valid = splits[[2]]
test = splits[[3]]
# dat_h2o = h2o.importFile(localH2O, path = ".../")

hyper_params <- list(
    activation=c("RectifierWithDropout","TanhWithDropout","MaxoutWithDropout"),
    hidden=list(  c(4000,4000) ),
    input_dropout_ratio=c( 0.2),
    l1 = c(1e-5, 1e-3, 0),
    l2 = c(1e-5, 1e-3, 0)
)

search_criteria = list(strategy = "RandomDiscrete", max_models = 50, seed=54321, stopping_rounds=5, stopping_tolerance=1e-2)

dl_random_grid <- h2o.grid(
    algorithm="deeplearning",
    grid_id = "dl_grid_random",
    training_frame=train,
    validation_frame=valid, 
    x=2:ncol(train), 
    y=1,
    epochs= 1,
    stopping_metric="RMSE",
    stopping_tolerance=1e-2,        ## stop when RMSE does not improve by >=1% for 5 scoring events
    stopping_rounds=5,            
    hyper_params = hyper_params,
    search_criteria = search_criteria
)          

grid <- h2o.getGrid("dl_grid_random" ,sort_by="RMSE",decreasing=FALSE)
grid
best_model <- h2o.getModel(grid@model_ids[[1]])


test_vals = h2o.predict(object = best_model, test)
plot(as.vector(as.numeric(test_vals)), as.vector(as.numeric(test$log2FoldChange) ))
abline(lm(as.matrix(test_vals) ~ as.matrix(test$log2FoldChange)), col="blue", lwd=2)
cor.test(as.vector(as.numeric(test_vals)), as.vector(as.numeric(test$log2FoldChange) ))


s <- proc.time()
H2O = h2o.deeplearning(x =2:ncol(train), #predictors
                       y = 1, #label 
                       training_frame = train,
                       activation = "RectifierWithDropout",
                       input_dropout_ratio = 0.2,
                       validation_frame = valid,
                       stopping_rounds = 5,
                       stopping_tolerance = 0.01,
                       stopping_metric = "RMSE",
                       seed = 54321,
                       l1 = 1e-5,
                       l2 = 1e-5,
                       hidden = c(3600 ,3600, 3600),
                       epochs = 500)

diff  = proc.time() - s

test_vals = h2o.predict(object = H2O, test)
plot(as.vector(as.numeric(test_vals)), as.vector(as.numeric(test$log2FoldChange) ))
abline(lm(as.matrix(test_vals) ~ as.matrix(test$log2FoldChange)), col="blue", lwd=2)
cor.test(as.vector(as.numeric(test_vals)), as.vector(as.numeric(test$log2FoldChange) ))


