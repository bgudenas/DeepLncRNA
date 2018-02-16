res = readRDS("./Data/Training_frames.rds")


library(h2o)
localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = 24, max_mem_size = "84G")
train_hex = as.h2o(res$train[ ,-1])
valid_hex = as.h2o(res$validate[ ,-1])
test_hex = as.h2o(res$validate[ ,-1])

hyper_params <- list(
  activation=c("RectifierWithDropout"),
  hidden =list( c(250,125,75), c(100,50,25 ), c(64, 32, 16), c(64,64,64) ),
                input_dropout_ratio=c(0, 0.1,0.2, 0.3, 0.4),
                l1 = c(1e-6,1e-5,1e-4, 1e-3, 0),
                l2 = c(1e-6,1e-5,1e-4, 1e-3, 0)
  )
  
  
  search_criteria = list(strategy = "RandomDiscrete", max_models = 800, seed=54321, stopping_rounds=10, stopping_tolerance=1e-3)
  
  dl_random_grid <- h2o.grid(
    algorithm="deeplearning",
    grid_id = "dl_grid",
    training_frame=train_hex,
    validation_frame = valid_hex,
    x=1:1583,
    y=1584,
    epochs= 400,
    stopping_metric="misclassification",
    stopping_tolerance=1e-3,        
    stopping_rounds=10,
    hyper_params = hyper_params,
    search_criteria = search_criteria
  )
  
  #print(dl_random_grid)
  grid <- h2o.getGrid("dl_grid",sort_by="max_per_class_error", decreasing=FALSE)
  grid2 <- h2o.getGrid("dl_grid",sort_by="accuracy", decreasing=TRUE)
  
  print(grid@summary_table[1:5,])
  print(grid2@summary_table[1:5,])
  best_model <- h2o.getModel(grid@model_ids[[1]])
  best_model
  
  #print(best_model@allparameters)
  print(h2o.performance(best_model, valid=T))
  
  print(h2o.performance(best_model, test_hex))
  h2o.saveModel(best_model, path="./Model/", force=TRUE)
  