res = readRDS("./Data/Training_frames.rds")
library(randomForest)
library(doParallel)
library(caret)
library(broom)

cluster <- makeCluster(detectCores()) 
registerDoParallel(cluster)

res$train$Loc = "Nuclear"
res$train$Loc[res$train$l2fc < 0] = "Cytosol"
res$train = res$train[res$train$l2fc < 0 | res$train$l2fc > 2.8 , ]
res$train$Loc = as.factor(res$train$Loc)
res$train = res$train[ ,-1]

res$validate$Loc = "Nuclear"
res$validate$Loc[res$validate$l2fc < 0] = "Cytosol"
res$validate = res$validate[res$validate$l2fc < 0 | res$validate$l2fc > 2.8 , ]
res$validate$Loc = as.factor(res$validate$Loc)

res$test$Loc = "Nuclear"
res$test$Loc[res$test$l2fc < 0] = "Cytosol"
res$test = res$test[res$test$l2fc < 0 | res$test$l2fc > 2.8 , ]
res$test$Loc = as.factor(res$test$Loc)


rf <- caret::train(Loc ~., data = res$train, method = "rf",
                           allowParralel = TRUE,
                           tuneLength = 5,
                           ntree = 128)


#rf_valid = sqrt(mean((valid_vals - res$validate$l2fc)^2)) 

#ptm <- proc.time()
svm_radial <- caret::train(Loc ~., data = res$train, method = "svmRadial",
                    preProcess = c("center", "scale"),
                    allowParralel = TRUE,
                    tuneLength = 5)
#proc.time() - ptm

#svm_vals = predict(svm_radial, res$validate)
#svm_valid = sqrt(mean((valid_vals - res$validate$l2fc)^2)) 

stopCluster(cluster)
registerDoSEQ()
 save.image("./Data/ML_models.RData")
 


library(pROC)
rf_preds = predict(rf, res$validate, type = "prob")
rf_test = predict(rf, res$validate)

rf.ROC <- pROC::roc(predictor=rf_preds$Cytosol,
               response=res$validate$Loc,
               levels=levels(res$validate$Loc) )



svm_preds = predict(svm_radial, res$validate,  type = "prob")

rf.ROC <- pROC::roc(predictor=rf_preds$Cytosol,
                    response=res$validate$Loc,
                    levels=levels(res$validate$Loc) )



library(h2o)
localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = 1, max_mem_size = "6G")
#DNN = h2o.loadModel("./Data/dl_grid_model_190")
DNN = h2o.loadModel("./Data/Models/dlgrid_model_151")
train_hex = as.h2o(res$train)
valid_hex = as.h2o(res$validate)
test_hex = as.h2o(res$test)


## Make plot of validation accuracies
rf_vals = predict(rf, res$validate)
rf_conf = confusionMatrix(rf_vals, res$validate$Loc)

svm_vals = predict(svm_radial, res$validate)
svm_conf = confusionMatrix(svm_vals, res$validate$Loc)

DNN_vals = predict(DNN, valid_hex)
DNN_conf = confusionMatrix(as.vector(DNN_vals$predict), res$validate$Loc)

GetMets = function(confMat){
  ## get acc, sens, spec from caret confusion matrix
  acc = confMat$overall[1]
  sens =  confMat$table[1,1]/(confMat$table[1,1] + confMat$table[2,1])
  spec = confMat$table[2,2]/(confMat$table[2,2] + confMat$table[1,2])
  mets = c(acc,sens,spec)
  names(mets) = c("Accuracy","Sensitivity","Specificity")
  return(mets)
}

valid_mets = data.frame(RF = GetMets(rf_conf), SVM = GetMets(svm_conf), DNN = GetMets(DNN_conf), Metrics = names(GetMets(DNN_conf))) %>% 
  tidyr::gather(valid_mets, Metrics)
colnames(valid_mets) = c("Metric","Model","Value")

ggplot(data = valid_mets, aes(as.factor(Model), Value, fill = Metric) ) +
  geom_bar( stat = "identity", position = "dodge") +
  scale_fill_brewer(palette = "Dark2") +
  coord_cartesian(ylim=c(0.65,0.8)) +
  xlab("Model") +
  theme_bw() + 
  theme(axis.text=element_text(size=13),
  axis.title=element_text(size=14,face="bold"), legend.title = element_text(size=12,face="bold"))

ggsave("./Figures/Model_Comparison.pdf")

svm_vals = predict(svm_radial, res$validate)
confusionMatrix(svm_vals, res$validate$Loc)

r
# caret::confusionMatrix(data = rf_test, res$test$Loc)
# 
# h2o.performance(DNN, test_hex)
