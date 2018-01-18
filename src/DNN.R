res = readRDS("./Data/Training_frames.rds")

res$train$Loc = "Nuclear"
res$train$Loc[res$train$l2fc <= -1] = "Cytosol"
res$train = res$train[res$train$l2fc <= -1 | res$train$l2fc >= 3 , ]
res$train$Loc = as.factor(res$train$Loc)

res$validate$Loc = "Nuclear"
res$validate$Loc[res$validate$l2fc <= -1] = "Cytosol"
res$validate = res$validate[res$validate$l2fc <= -1 | res$validate$l2fc >= 3 , ]
res$validate$Loc = as.factor(res$validate$Loc)

library(h2o)
localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = 24, max_mem_size = "64G")
train_hex = as.h2o(res$train[ ,-1])
valid_hex = as.h2o(res$validate[ ,-1])


hyper_params <- list(
  activation=c("RectifierWithDropout"),
  hidden =list( c(2000,2000,2000 ), c(1500,1500,1500 ), c(1000,1000,1000), c(500,500,500) ),
                input_dropout_ratio=c(0, 0.1,0.2),
                l1 = c(1e-6,1e-3, 1e-4, 0),
                l2 = c(1e-6,1e-3, 1e-4, 0)
  )
  
  
  search_criteria = list(strategy = "RandomDiscrete", max_models = 300, seed=54321, stopping_rounds=5, stopping_tolerance=1e-2)
  
  dl_random_grid <- h2o.grid(
    algorithm="deeplearning",
    grid_id = "dl_grid",
    training_frame=train_hex,
    validation_frame = valid_hex,
    x=1:1557,
    y=1558,
    epochs= 5,
    stopping_metric="AUC",
    stopping_tolerance=1e-2,        ## stop when RMSE does not improve by >=1% for 5 scoring events
    stopping_rounds=5,
    hyper_params = hyper_params,
    search_criteria = search_criteria
  )
  
  #print(dl_random_grid)
  grid <- h2o.getGrid("dl_grid",sort_by="AUC",decreasing=FALSE)
  
  print(grid@summary_table[1:5,])
  best_model <- h2o.getModel(grid@model_ids[[1]])
  best_model
  
  #print(best_model@allparameters)
  print(h2o.performance(best_model, valid=T))
  
  
  h2o.saveModel(best_model, path="./Model/DNN", force=TRUE)
  