library(h2o)
set.seed(54321)
localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = -1, max_mem_size = "15G")
# Import data -------------------------------------------------------------
df = read.csv("./Data/lncRNA_featureset.csv", row.names = 1)
filt = -c(1,2,4,5,6,7,8,13)
df = df[ , filt]
dat_h2o = as.h2o(df)


splits = h2o.splitFrame(dat_h2o, ratios =  c(0.7), seed = 54321)
train_hex = splits[[1]]
test_hex = splits[[2]]


train <- train_hex[,-1]
test  <- test_hex [,-1]


    get_stacked_ae_array <- function(training_data,layers,args){  
        vector <- c()
        index = 0
        for(i in 1:length(layers)){    
            index = index + 1
            ae_model <- do.call(h2o.deeplearning, 
                                modifyList(list(x=names(training_data),
                                                training_frame=training_data,
                                                autoencoder=T,
                                                hidden=layers[i]),
                                           args))
            training_data = h2o.deepfeatures(ae_model,training_data,layer=1)
            
            names(training_data) <- gsub("DF", paste0("L",index,sep=""), names(training_data)) 
            vector <- c(vector, ae_model)    
        }
        vector
    }
    
    # this function returns final encoded contents
    apply_stacked_ae_array <- function(data,ae){
        index = 0
        for(i in 1:length(ae)){
            index = index + 1
            data = h2o.deepfeatures(ae[[i]],data,layer=1)
            names(data) <- gsub("DF", paste0("L",index,sep=""), names(data)) 
        }
        data
    }

    ## Build reference model on full dataset and evaluate it on the test set
    # model_ref <- h2o.deeplearning(training_frame=train_hex, x=2:(ncol(train_hex)), y=response, hidden=c(1000,500,250,70), epochs=1)
    # p_ref <- h2o.performance(model_ref, test_hex)
    # p_ref
    
    ## Now build a stacked autoencoder model with three stacked layer AE models
    ## First AE model will compress the 717 non-const predictors into 200
    ## Second AE model will compress 200 into 100
    ## Third AE model will compress 100 into 50
    layers <- c(2000, 1000, 500)
    args <- list(activation= "Tanh", epochs=100, input_dropout_ratio = 0.1)
    ae <- get_stacked_ae_array(train, layers, args)
    ## Now compress the training/testing data with this 3-stage set of AE models
    train_compressed <- apply_stacked_ae_array(train, ae)
    test_compressed <- apply_stacked_ae_array(test, ae)
    
    ## Build a simple model using these new features (compressed training data) and evaluate it on the compressed test set.
    train_w_resp <- h2o.cbind(train_compressed, train_hex[,1])
    test_w_resp <- h2o.cbind(test_compressed, test_hex[,1])
    save.image("SDA.RData")  
    
    model_on_compressed_data <- h2o.deeplearning(training_frame = train_w_resp,
                                                 x=1:(ncol(train_w_resp)-1),
                                                 y=ncol(train_w_resp),
                                                 validation_frame = test_w_resp,
                                                 stopping_metric = "RMSE",
                                                 stopping_tolerance = 0.001,
                                                 stopping_rounds = 10,
                                                 hidden=c(200,200, 200),
                                                 activation = "RectifierWithDropout",
                                                 input_dropout_ratio = 0.2,
                                                 l2=1e-5,
                                                 model_id = "Sup_DNN",
                                                 epochs=1000)
    
     h2o.performance(model_on_compressed_data, test_w_resp)
     h2o.saveModel(object = model_on_compressed_data, path = getwd())
     
    train_vals = h2o.predict(object = model_on_compressed_data, train_w_resp)
    test_vals = h2o.predict(object = model_on_compressed_data, test_w_resp)
    plot(as.vector(as.numeric(test_vals)), as.vector(as.numeric(test_w_resp[ ,ncol(test_w_resp)]) ) )
    abline(lm(as.matrix(test_vals) ~ as.matrix(test$log2FoldChange)), col="blue", lwd=2)
    cor.test(as.vector(as.numeric(test_vals)), as.vector(as.numeric(test_w_resp[ ,ncol(test_w_resp)]) ))
    
    trainRF = as.data.frame(train_w_resp)
    testRF = as.data.frame(test_w_resp)
    #AUC:  0.7156318 : .89 cytosol error
    library(randomForest)
    rf = randomForest(log2FoldChange ~ .,
                      data = trainRF,
                      ntree=2000)
    
    test_valsRF = predict(rf, testRF)
    
    cor.test(as.vector(as.numeric(test_valsRF)), as.vector(as.numeric(test_w_resp[ ,ncol(test_w_resp)]) ))
    RMSE <- sqrt(mean((test_valsRF - as.vector(as.numeric(test_w_resp[ ,ncol(test_w_resp)]) )) ^2))
    
    plot(test_vals, as.vector(as.numeric(test_w_resp[ ,ncol(test_w_resp)]) ), ylim = range(-2, 4), xlim=range(-1, 3) )
    abline(lm(as.vector(as.numeric(test_w_resp[ ,ncol(test_w_resp)]) ) ~ test_vals), col="blue", lwd=2)
    
  
    save.image("./Data/SDA_models.RData")
    
    
# build training frame ----------------------------------------------------
    train_vals = h2o.predict(object = model_on_compressed_data, train_w_resp)
    train_RF = predict(rf, trainRF)
    
    trains = as.data.frame(matrix(nrow= nrow(train_vals), data = NA))
    trains$DNN = as.vector(train_vals)
    trains$SDA_RF = as.vector(train_RF)
    trains$L2FC = as.vector(train_w_resp[ ,ncol(train_w_resp)])
    write.csv(trains, "./Data/Training_preds.csv")
    
# build testing frame -----------------------------------------------------
    test_vals = h2o.predict(object = model_on_compressed_data, test_w_resp)
    test_valsRF = predict(rf, testRF)
    
    tests = as.data.frame(matrix(nrow= nrow(test_vals), data = NA))
    tests$DNN = as.vector(test_vals)
    tests$SDA_RF = as.vector(test_valsRF)
    tests$L2FC = as.vector(test_w_resp[ ,ncol(test_w_resp)])
    write.csv(tests, "./Data/Tests_preds.csv")

