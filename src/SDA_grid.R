source("./src/SDA_funcs.R")
library(h2o)
set.seed(54321)

# Import data -------------------------------------------------------------
df = readRDS("./Data/Training_frames.rds")


localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = -1, max_mem_size = "15G")


train_hex = as.h2o(df$train[ ,-1])
#train <- train_hex[,-1]

# Build SDAs --------------------------------------------------------------
# test using diff architectures
layer1 = c(1500, 750, 375)
layer_ov =c(2000,1000,500)
layer_slim =c(750,350,150)

args1 <- list(activation= "Tanh", epochs=100, input_dropout_ratio = 0.1, mini_batch_size = 1)
args2 <- list(activation= "Tanh", epochs=100, input_dropout_ratio = 0.25, mini_batch_size = 1)
args3 <- list(activation= "Tanh", epochs=100, input_dropout_ratio = 0.4, mini_batch_size = 1)

#layer1
ael1 <- get_stacked_ae_array(train_hex, layer1, args1)
ael2 <- get_stacked_ae_array(train_hex, layer1, args2)
ael3 <- get_stacked_ae_array(train_hex, layer1, args3)

train_ael1 <- apply_stacked_ae_array(train_hex, ael1)
train_ael2 <- apply_stacked_ae_array(train_hex, ael2)
train_ael3 <- apply_stacked_ae_array(train_hex, ael3)

#layerov
aeov1 <- get_stacked_ae_array(train_hex, layer_ov, args1)
aeov2 <- get_stacked_ae_array(train_hex, layer_ov, args2)
aeov3 <- get_stacked_ae_array(train_hex, layer_ov, args3)

train_ov1 <- apply_stacked_ae_array(train_hex, aeov1)
train_ov2 <- apply_stacked_ae_array(train_hex, aeov2)
train_ov3 <- apply_stacked_ae_array(train_hex, aeov3)

#layer_slim
aesl1 <- get_stacked_ae_array(train_hex, layer_slim, args1)
aesl2 <- get_stacked_ae_array(train_hex, layer_slim, args2)
aesl3 <- get_stacked_ae_array(train_hex, layer_slim, args3)

train_sl1 <- apply_stacked_ae_array(train_hex, aesl1)
train_sl2 <- apply_stacked_ae_array(train_hex, aesl2)
train_sl3 <- apply_stacked_ae_array(train_hex, aesl3)

#save.image("./SDA_grid.RData")

sink('./Data/SDA_gridsearch.txt')
SDA_rf(train_ael1)
SDA_rf(train_ael2)
SDA_rf(train_ael3)

SDA_rf(train_ov1)
SDA_rf(train_ov2)
SDA_rf(train_ov3)

SDA_rf(train_sl1)
SDA_rf(train_sl2)
SDA_rf(train_sl3)
sink()
h2o.shutdown()
