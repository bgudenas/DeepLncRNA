source("./src/SDA_funcs.R")
library(h2o)
set.seed(54321)

# Import data -------------------------------------------------------------
df = readRDS("../Data/Unsupervised_train.rds")
Super = readRDS("../Data/Supervised_train.rds")
df$transcript_biotype = as.factor(df$transcript_biotype)
Super$transcript_biotype = as.factor(Super$transcript_biotype)

df = df[ , -c(3, 5)]
Super = Super[ , -5]


localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = 2, max_mem_size = "6G")

dat_h2o = as.h2o(df)
dat_super = as.h2o(Super)
splits = h2o.splitFrame(dat_super, ratios =  c(0.7), seed = 54321)

train_hex = splits[[1]]
test_hex = splits[[2]]

train <- train_hex[,-1]
test  <- test_hex [,-1]

# Build SDAs --------------------------------------------------------------
# test using diff architectures
layer1 = c(1500,700,300)
layer_ov =c(2000,1000,500)
layer_slim =c(900,300,100)

args1 <- list(activation= "Tanh", epochs=1, input_dropout_ratio = 0.1)
args2 <- list(activation= "Tanh", epochs=1, input_dropout_ratio = 0.25)

#layer1
ael1 <- get_stacked_ae_array(dat_h2o, layer1, args1)
ael2 <- get_stacked_ae_array(dat_h2o, layer1, args2)

train_ael1 <- apply_stacked_ae_array(train, ael1)
train_ael2 <- apply_stacked_ae_array(train, ael2)

#layerov
aeov1 <- get_stacked_ae_array(dat_h2o, layer_ov, args1)
aeov2 <- get_stacked_ae_array(dat_h2o, layer_ov, args2)

train_ov1 <- apply_stacked_ae_array(train, aeov1)
train_ov2 <- apply_stacked_ae_array(train, aeov2)

#layer_slim
aesl1 <- get_stacked_ae_array(dat_h2o, layer_slim, args1)
aesl2 <- get_stacked_ae_array(dat_h2o, layer_slim, args2)

train_sl1 <- apply_stacked_ae_array(train, aesl1)
train_sl2 <- apply_stacked_ae_array(train, aesl2)

save.image("../Data/SDA_grid.RData")


SDA_rf(train_ael1)
SDA_rf(train_ael2)

SDA_rf(train_ov1)
SDA_rf(train_ov2)

SDA_rf(train_sl1)
SDA_rf(train_sl1)
