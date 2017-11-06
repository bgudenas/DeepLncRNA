
# AE ----------------------------------------------------------------------
df = readRDS("./Data/Supervised_train.rds")
df$transcript_biotype = as.factor(df$transcript_biotype)
df = df[ ,-5]
Loc = rep("", nrow(df))
Loc[df$b < -1] = "Cytosol"
Loc[df$b > 3] = "Nuclear"
df = cbind(Loc, df)

set.seed(54321)
inTrain = sample(1:nrow(df), .7*nrow(df) )
train = df[inTrain, ]
test = df[-inTrain, ]

train = train[train$Loc != "", ]
train$Loc = droplevels(train$Loc)

test = test[test$Loc != "", ]
test$Loc = droplevels(test$Loc)


library(h2o)
localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = -1, max_mem_size = "15G")
train_hex = as.h2o(train)
test_hex = as.h2o(test)

hid1 = c(1000,500,200,500,1000)
hid = c(900,400,200,400,900)
hid2 =c(1000,700,500,300)
AE = h2o.deeplearning(
  x = 3:1557,
  training_frame = train_hex,
  activation="Rectifier",
  max_w2 = 10,
  hidden = hid1,
  #input_dropout_ratio = 0.1,
  ignore_const_cols=F,
  epochs = 50,
  autoencoder = TRUE,
  model_id="AE")

# train_sup = train[train$Loc != "", ]
# train_sup$Loc = droplevels(train_sup$Loc)
# test_sup = test[test$Loc != "", ]
# test_sup$Loc = droplevels(test_sup$Loc)
# 
# train_sup = as.h2o(train_sup)
# test_sup = as.h2o(test_sup)

pretrained_model <- h2o.deeplearning(x=3:1557, 
                                     y=1, 
                                     model_id = "DNN",
                                     activation = "RectifierWithDropout",
                                     training_frame=train_hex,
                                     validation_frame=test_hex,
                                     hidden= hid1,
                                     epochs=200, 
                                    #balance_classes = TRUE,
                                    input_dropout_ratio = 0.1,
                                    l1 = 1e-4,
                                    l2 = 1e-3 ,
                                    seed = 54321) #,
                                    #pretrained_autoencoder="AE")
pretrained_model

h2o.saveModel(pretrained_model, path = paste0(getwd(),"/Data"), force = TRUE)

hyper_params <- list(
  activation=c("RectifierWithDropout", "Rectifier"),
  hidden =list(c(1000,500,250,500,1000), c(2000,1000,500,250), c(1000,750,500,250), c(2000,1000,500,250,125) ),
  input_dropout_ratio=c( 0.1,0.2, 0.3, 0.4),
  l1 = c(1e-6,1e-5, 1e-4, 0),
  l2 = c(1e-6,1e-5, 1e-4, 0)
)


search_criteria = list(strategy = "RandomDiscrete", max_models = 200, seed=54321, stopping_rounds=5, stopping_tolerance=1e-2)

dl_random_grid <- h2o.grid(
  algorithm="deeplearning",
  grid_id = "dl_grid_random",
  training_frame=train_hex,
  nfolds = 3,
  x=3:1557, 
  y=1,
  epochs= 1,
  stopping_metric="AUC",
  stopping_tolerance=1e-2,        ## stop when RMSE does not improve by >=1% for 5 scoring events
  stopping_rounds=5,            
  hyper_params = hyper_params,
  search_criteria = search_criteria
)          

print(dl_random_grid)

  
  