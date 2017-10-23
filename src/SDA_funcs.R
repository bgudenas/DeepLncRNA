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
    as.data.frame(data)
}


SDA_rf = function(trains){
    print(dim(trains))
    resp_num = ncol(trains)+1
    trains = as.h2o(trains)
    trains = h2o.cbind(trains, train_hex$b)
    
    rf1 = h2o.randomForest(training_frame = trains,
                           ntrees = 100,
                           y = resp_num,
                           nfolds = 5)
    return(print(rf1))
}
    
    