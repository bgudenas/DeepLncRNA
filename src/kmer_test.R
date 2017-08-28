
# Test K for K-mers -------------------------------------------------------

library(h2o)
set.seed(54321)

# Import data -------------------------------------------------------------
df = read.csv("./Data/lncRNA_featureset.csv", row.names = 1)
filt = -c(1,2,4,5,6,7,8,13)
df = df[ , filt]


for (i in 1:nrow(df)){
  T_length = df$transcript_length[i]
  vec = as.vector(df[i, 6:5381 ]/T_length)
  df[i, 6:5381 ] = vec
}

log_con = file("kmer_test.txt")

# Start H2O cluster -------------------------------------------------------
localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = -1, max_mem_size = "14G")

dat_h2o = as.h2o(df)
splits = h2o.splitFrame(dat_h2o, ratios =  c(0.6, 0.2), seed = 54321)
train = splits[[1]]
valid = splits[[2]]
test = splits[[3]]



sixs = 6:(5+4096)
fives = (4102):(4101 + 1024)
fours = 5126:(5125+256)



H2O_6 = h2o.deeplearning(x = sixs, #predictors
                       y = 1, #label 
                       training_frame = train,
                       activation = "RectifierWithDropout",
                       input_dropout_ratio = 0.2,
                       validation_frame = valid,
                       stopping_rounds = 5,
                       stopping_metric = "RMSE",
                       stopping_tolerance = 0.01,
                       seed = 54321,
                       l2 = 1e-5,
                       hidden = c(4096, 2048, 1024),
                       epochs = 300)
cat("6mers------------------------------",  file = log_con)
cat(print(H2O_6), file = log_con)

test_vals = h2o.predict(object = H2O_6, test)
cor = cor.test(as.vector(as.numeric(test_vals)), as.vector(as.numeric(test$log2FoldChange) ))

p = h2o.performance(H2O_6, test)
cat(p, file = log_con)
cat(cor, file = log_con)

# 5mers -------------------------------------------------------------------
H2O_5 = h2o.deeplearning(x = fives, #predictors
                         y = 1, #label 
                         training_frame = train,
                         activation = "RectifierWithDropout",
                         input_dropout_ratio = 0.2,
                         validation_frame = valid,
                         stopping_rounds = 5,
                         stopping_metric = "RMSE",
                         stopping_tolerance = 0.01,
                         seed = 54321,
                         l2 = 1e-5,
                         hidden = c(1024, 512, 256),
                         epochs = 300)
cat("5mers------------------------------", file = log_con)
cat(print(H2O_5), file = log_con)

test_vals = h2o.predict(object = H2O_5, test)
cor = cor.test(as.vector(as.numeric(test_vals)), as.vector(as.numeric(test$log2FoldChange) ))

p = h2o.performance(H2O_5, test)
cat(p, file = log_con)
cat(cor, file = log_con)

# 4mers -------------------------------------------------------------------
H2O_4 = h2o.deeplearning(x = fours, #predictors
                         y = 1, #label 
                         training_frame = train,
                         activation = "RectifierWithDropout",
                         input_dropout_ratio = 0.2,
                         validation_frame = valid,
                         stopping_rounds = 5,
                         stopping_metric = "RMSE",
                         stopping_tolerance = 0.01,
                         seed = 54321,
                         l2 = 1e-5,
                         hidden = c(256, 128, 64),
                         epochs = 300)
cat("4mers------------------------------", file = log_con)
cat(print(H2O_4), file = log_con)

test_vals = h2o.predict(object = H2O_4, test)
cor = cor.test(as.vector(as.numeric(test_vals)), as.vector(as.numeric(test$log2FoldChange) ))

p = h2o.performance(H2O_4, test)
cat(p, file = log_con)
cat(cor, file = log_con)

