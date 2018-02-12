res = readRDS("./Data/Training_frames.rds")
library(randomForest)
library(doParallel)
library(caret)

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

svm_radial <- caret::train(Loc ~.,
                           data = res$train,
                           method = "svmRadial",
                           preProcess = c("center", "scale"),
                           trControl = trainControl(classProbs = TRUE),
                           allowParralel = TRUE,
                           tuneLength = 3)

rf <- caret::train(Loc ~ ., data = res$train, method = "rf",
                           allowParralel = TRUE,
                           tuneLength = 1,
                           ntree = 101)



stopCluster(cluster)
registerDoSEQ()
 save.image("./Data/ML_models.RData")
 load(file = "./Data/ML_models.RData")
 


library(pROC)
library(h2o)
 library(caret)
 library(dplyr)
localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = 1, max_mem_size = "6G")
#DNN = h2o.loadModel("./Data/dl_grid_model_190")
DNN = h2o.loadModel("./Data/Models/dlgrid_model_483")
train_hex = as.h2o(res$train)
valid_hex = as.h2o(res$validate)
test_hex = as.h2o(res$test)


## Make plot of validation accuracies
rf_vals = predict(rf, res$validate)
rf_conf = confusionMatrix(rf_vals, res$validate$Loc, positive = "Nuclear")

svm_vals = predict(svm_radial, res$validate)
svm_conf = confusionMatrix(svm_vals, res$validate$Loc, positive = "Nuclear")

DNN_vals = predict(DNN, valid_hex)
DNN_conf = confusionMatrix(as.vector(DNN_vals$predict), res$validate$Loc, positive = "Nuclear")


GetMets = function(confMat){
  ## get acc, sens, spec from caret confusion matrix
  acc = confMat$overall[1]
  spec =  confMat$table[1,1]/(confMat$table[1,1] + confMat$table[2,1])
  sens = confMat$table[2,2]/(confMat$table[2,2] + confMat$table[1,2])
  mets = c(acc,sens,spec)
  names(mets) = c("Accuracy","Sensitivity","Specificity")
  return(mets)
}

valid_mets = data.frame(RF = GetMets(rf_conf), SVM = GetMets(svm_conf), DeepLocal = GetMets(DNN_conf), Metrics = names(GetMets(DNN_conf))) %>% 
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



# caret::confusionMatrix(data = rf_test, res$test$Loc)
# make ROC plot on test set

dnn_test = predict(DNN, test_hex)
rf_test = predict(rf, res$test, "prob")
svm_test = predict(svm_radial, res$test, "prob")
pdf("./Figures/AUROC.pdf")


dnn_roc = pROC::roc(predictor=as.vector(dnn_test$Nuclear),
                    response=res$test$Loc,
                    levels=levels(res$test$Loc) )

plot(dnn_roc, col = "blue", cex.lab = 1.3, cex.axis = 1, lwd=3)

rf_roc =pROC::roc(predictor=as.vector(rf_test$Nuclear),
                  response=res$test$Loc,
                  levels=levels(res$test$Loc) )

plot(rf_roc, add = TRUE, col = "darkorange",lty=2, lwd = 3)

svm_roc = pROC::roc(predictor=as.vector(svm_test$Nuclear),
                     response=res$test$Loc,
                     levels=levels(res$test$Loc) )
plot(svm_roc, add = TRUE, col = "magenta",lty =3, lwd = 3)


legend(0.35, 0.25, legend=c("DNN", "RF", "SVM","Random Guess"),
       col=c("blue", "darkorange","magenta","grey"), lty = c(1, 2,3,1), lwd = c(3,3,3,2), cex = 1.2)
dev.off()

rf_test = predict(rf, res$test)
rf_conf = confusionMatrix(rf_test, res$test$Loc, positive = "Nuclear")

svm_test = predict(svm_radial, res$test)
svm_conf = confusionMatrix(svm_test, res$test$Loc, positive = "Nuclear")

DNN_test = predict(DNN, test_hex)
DNN_conf = confusionMatrix(as.vector(DNN_test$predict), res$test$Loc, positive = "Nuclear")

test_mets = data.frame(t(data.frame(RF = GetMets(rf_conf), SVM = GetMets(svm_conf), DNN = GetMets(DNN_conf)) ))
test_mets$AUC = c(rf_roc$auc, svm_roc$auc, dnn_roc$auc)

rf_mcc = mltools::mcc(preds = as.numeric(rf_test)-1, actuals = as.numeric(res$test$Loc)-1)
svm_mcc = mltools::mcc(preds = as.numeric(svm_test)-1, actuals = as.numeric(res$test$Loc)-1)
dnn_mcc = mltools::mcc(preds = as.numeric(as.factor(as.vector(DNN_test$predict)))-1, actuals = as.numeric(res$test$Loc)-1)

test_mets$MCC = c(rf_mcc, svm_mcc, dnn_mcc)
write.csv(round(test_mets, 3), "./Data/test_metrics.csv")
#   
# test_mets
# RF       SVM       DNN     Metrics
# Accuracy    0.7162162 0.7132132 0.7327327    Accuracy
# Sensitivity 0.7071742 0.7057101 0.6617862 Sensitivity
# Specificity 0.7257319 0.7211094 0.8073960 Specificity

h2o.shutdown(FALSE)
