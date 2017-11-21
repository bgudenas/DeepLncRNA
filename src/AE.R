library(h2o)
set.seed(54321)

# Import data -------------------------------------------------------------
df = readRDS("../Data/Unsupervised_train.rds")
Super = readRDS("../Data/Supervised_train.rds")
colnames(Super)[2] = "transcript_biotype"
df = df[ , -c(3, 5]
Super = Super[ , -5]


# Loc = cbind(local = rep("", nrow(df)), df)
# Loc$local = as.character(Loc$local)
# Loc$local[Loc$log2FoldChange <= -0.5] = "Cytosol"
# Loc$local[Loc$log2FoldChange >= 2] = "Nuclear"
# Loc = Loc[ ,-2]
# Loc$local = as.factor(Loc$local)
# 
# Loc_un = which(Loc$local == "")
# Start H2O cluster -------------------------------------------------------

localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = 16, max_mem_size = "50G")

dat_h2o = as.h2o(df)
dat_super = as.h2o(Super)
splits = h2o.splitFrame(dat_super, ratios =  c(0.7), seed = 54321)

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
    layers <- c(2500, 1000, 500)
    args <- list(activation= "Tanh", epochs=1000, input_dropout_ratio = 0.2)
    ae <- get_stacked_ae_array(dat_h2o, layers, args)
    #model_path = h2o.saveModel(object = ae, path = getwd(), force = TRUE)
    
    ## Now compress the training/testing data with this 3-stage set of AE models
    train_compressed <- apply_stacked_ae_array(train, ae)
    test_compressed <- apply_stacked_ae_array(test, ae)
#     
#     ## Build a simple model using these new features (compressed training data) and evaluate it on the compressed test set.
#     train_w_resp <- h2o.cbind(train_compressed, train_hex[,1])
#     test_w_resp <- h2o.cbind(test_compressed, test_hex[,1])
#     
#     model_on_compressed_data <- h2o.deeplearning(training_frame = train_w_resp,
#                                                  x=1:(ncol(train_w_resp)-1),
#                                                  y=ncol(train_w_resp),
#                                                  validation_frame = test_w_resp,
#                                                  stopping_metric = "RMSE",
#                                                  stopping_tolerance = 0.01,
#                                                  hidden=c(500),
#                                                  activation = "Rectifier",
#                                                  input_dropout_ratio = 0.2,
#                                                  # balance_classes = TRUE,
#                                                  epochs=1000,
# 						model_id = "SUP_AE")
save.image("../Data/SDA.RData")
# h2o.saveModel("SUP_AE",path = getwd(),  force = TRUE) 
