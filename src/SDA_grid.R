source("./SDA_funcs.R")
library(h2o)
set.seed(54321)

# Import data -------------------------------------------------------------
df = readRDS("./Unsupervised_train.rds")
Super = readRDS("./Training_set.rds")

localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = 20, max_mem_size = "78G")

unsup = as.h2o(df[ ,-1])

train_hex = as.h2o(Super)
train <- train_hex[,-1]

# Build SDAs --------------------------------------------------------------
# test using diff architectures
layer1 = c(1500,1500,1500,1500)
layer_ov =c(4000,2000,1000,500)
layer_slim =c(3000,1500,750,350)

args1 <- list(activation= "Tanh", epochs=10, input_dropout_ratio = 0.1, mini_batch_size = 10)
args2 <- list(activation= "Tanh", epochs=10, input_dropout_ratio = 0.25, mini_batch_size = 10)
args3 <- list(activation= "Tanh", epochs=10, input_dropout_ratio = 0.4, mini_batch_size = 10)

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

save.image("./SDA_grid.RData")


SDA_rf(train_ael1, args1)
SDA_rf(train_ael2, args2)
SDA_rf(train_ael3, args3)

SDA_rf(train_ov1, args1)
SDA_rf(train_ov2, args2)
SDA_rf(train_ov3, args3)

SDA_rf(train_sl1, args1)
SDA_rf(train_sl2, args2)
SDA_rf(train_sl3, args3)
