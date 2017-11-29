source("./SDA_funcs.R")
library(h2o)
set.seed(54321)

# Import data -------------------------------------------------------------
df = readRDS("../Data/Unsupervised_train.rds")
Super = readRDS("../Data/Supervised_train.rds")
df$transcript_biotype = as.factor(df$transcript_biotype)
Super$transcript_biotype = as.factor(Super$transcript_biotype)

df = df[ , -c(3, 5)]
Super = Super[ , -5]


localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = 40, max_mem_size = "148G")

#dat_h2o = as.h2o(df)

inTrain = sample(1:nrow(Super), .7*nrow(Super) )

train_raw = Super[inTrain, ]
test_raw = Super[-inTrain, ]

## Remove any testing instances from unsupervised DF for SDA
df = df[is.na(match(rownames(df), rownames(test_raw))), ]
unsup = as.h2o(df)

train_hex = as.h2o(train_raw)
test_hex = as.h2o(test_raw)


train <- train_hex[,-1]
test  <- test_hex [,-1]

# Build SDAs --------------------------------------------------------------
# test using diff architectures
layer1 = c(3000,1500,750,350)
layer_ov =c(2000,1000,500,250)
layer_slim =c(500,500,500,500)

args1 <- list(activation= "Tanh", epochs=10, input_dropout_ratio = 0)
args2 <- list(activation= "Tanh", epochs=10, input_dropout_ratio = 0.1)
args3 <- list(activation= "Tanh", epochs=10, input_dropout_ratio = 0.25)

#layer1
ael1 <- get_stacked_ae_array(unsup, layer1, args1)
ael2 <- get_stacked_ae_array(unsup, layer1, args2)
ael3 <- get_stacked_ae_array(unsup, layer1, args3)

train_ael1 <- apply_stacked_ae_array(train, ael1)
train_ael2 <- apply_stacked_ae_array(train, ael2)
train_ael3 <- apply_stacked_ae_array(train, ael3)

#layerov
aeov1 <- get_stacked_ae_array(unsup, layer_ov, args1)
aeov2 <- get_stacked_ae_array(unsup, layer_ov, args2)
aeov3 <- get_stacked_ae_array(unsup, layer_ov, args3)

train_ov1 <- apply_stacked_ae_array(train, aeov1)
train_ov2 <- apply_stacked_ae_array(train, aeov2)
train_ov3 <- apply_stacked_ae_array(train, aeov3)

#layer_slim
aesl1 <- get_stacked_ae_array(unsup, layer_slim, args1)
aesl2 <- get_stacked_ae_array(unsup, layer_slim, args2)
aesl3 <- get_stacked_ae_array(unsup, layer_slim, args3)

train_sl1 <- apply_stacked_ae_array(train, aesl1)
train_sl2 <- apply_stacked_ae_array(train, aesl2)
train_sl3 <- apply_stacked_ae_array(train, aesl3)

save.image("../Data/SDA_grid.RData")


SDA_rf(train_ael1, args1)
SDA_rf(train_ael2, args2)
SDA_rf(train_ael3, args3)

SDA_rf(train_ov1, args1)
SDA_rf(train_ov2, args2)
SDA_rf(train_ov3, args3)

SDA_rf(train_sl1, args1)
SDA_rf(train_sl2, args2)
SDA_rf(train_sl3, args3)
