
# RF model on SDA features ------------------------------------------------
library(randomForest)

load(file="./Data/SDA_best.RData")
Loc = rep("", nrow(Super))
Loc[Super$b < -1] = "Cytosol"
Loc[Super$b > 3] = "Nuclear"
Super$Loc = as.factor(Loc)

set.seed(54321)
inTrain = sample(1:nrow(Super), .7*nrow(Super) )
train = Super[inTrain, ]
test = Super[-inTrain, ]


train_compressed$Loc = train$Loc
train_compressed = train_compressed[train_compressed$Loc != "", ]
train_compressed$Loc = droplevels(train_compressed$Loc)

test_compressed$Loc = test$Loc
test_compressed = test_compressed[test_compressed$Loc != "", ]
test_compressed$Loc = droplevels(test_compressed$Loc)

# library(randomForest)
# rf = randomForest(x = train_compressed[ , -151],
#                   y = train_compressed$Loc,
#                   ntree = 100)

library(caret)
tunegrid <- expand.grid(.mtry=seq(8,20,by = 2))
rf <- train(Loc ~ ., data = train_compressed, 
            method = "rf",
            tuneGrid = tunegrid,
            ntree = 100)

# mtry  Accuracy   Kappa    
# 8    0.7321926  0.4194834
# 10    0.7315515  0.4187037
# 12    0.7298732  0.4159250
# 14    0.7359974  0.4284726
# 16    0.7288430  0.4141539
# 18    0.7304584  0.4165074
# 20    0.7300365  0.4160999


rf = randomForest(x = train_compressed[ , -151],
                                     y = train_compressed$Loc,
                                     mtry = 14,
                                     ntree = 500)
saveRDS(rf, "./Data/SDA_RF.rds")
train_vals = predict(rf, train_compressed, predict.all = TRUE)
test_vals = predict(rf, test_compressed, predict.all = TRUE)
rf_preds = list(train_vals, test_vals)
saveRDS(rf_preds, "./Data/SDA_RF_preds.rds")

confusionMatrix(data = test_vals$aggregate, reference = test_compressed$Loc, mode = "everything")
# 
# Confusion Matrix and Statistics
# 
# Reference
# Prediction Cytosol Nuclear
# Cytosol     176      70
# Nuclear     141     418
# 
# Accuracy : 0.7379         
# 95% CI : (0.7061, 0.768)
# No Information Rate : 0.6062         
# P-Value [Acc > NIR] : 2.727e-15      
# 
# Kappa : 0.4286         
# Mcnemar's Test P-Value : 1.443e-06      
#                                          
#             Sensitivity : 0.5552         
#             Specificity : 0.8566         
#          Pos Pred Value : 0.7154         
#          Neg Pred Value : 0.7478         
#               Precision : 0.7154         
#                  Recall : 0.5552         
#                      F1 : 0.6252         
#              Prevalence : 0.3938         
#          Detection Rate : 0.2186         
#    Detection Prevalence : 0.3056         
#       Balanced Accuracy : 0.7059         
#                                          
#        'Positive' Class : Cytosol      
