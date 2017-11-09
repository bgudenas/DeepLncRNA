
# SDA Models --------------------------------------------------------------
load(file="./Data/SDA_best.RData")
#Data
set.seed(54321)
inTrain = sample(1:nrow(Super), .7*nrow(Super) )
train = Super[inTrain, ]
test = Super[-inTrain, ]
train_compressed = cbind(train$b, train_compressed)
test_compressed = cbind(test$b, test_compressed)


localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = -1, max_mem_size = "15G")
train_hex = as.h2o(train_compressed)
test_hex = as.h2o(test_compressed)
colnames(test_hex)[1] = "b"
colnames(train_hex)[1] = "b"

hyper_params <- list(
  activation=c("RectifierWithDropout", "Rectifier"),
  hidden =list( c(500,500,500), c(1000,500,250), c(300,100,300), c(250,250,250) ),
  input_dropout_ratio=c( 0.1,0.2, 0.4),
  l1 = c(1e-6,1e-5, 1e-4, 0),
  l2 = c(1e-6,1e-5, 1e-4, 0)
)

search_criteria = list(strategy = "RandomDiscrete", max_models = 200, seed=54321, stopping_rounds=5, stopping_tolerance=1e-2)

dl_random_grid <- h2o.grid(
  algorithm="deeplearning",
  grid_id = "dl_grid_random",
  training_frame=train_hex,
  nfolds = 3,
  x=2:1557, 
  y=1,
  epochs= 1,
  stopping_metric="deviance",
  stopping_tolerance=1e-2,        ## stop when deviance does not improve by >=1% for 5 scoring events
  stopping_rounds=5,            
  hyper_params = hyper_params,
  search_criteria = search_criteria
)          
