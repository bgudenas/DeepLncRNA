
# Test ML algos on raw data -----------------------------------------------
df = readRDS("./Data/Supervised_train.rds")
df$transcript_biotype = as.factor(df$transcript_biotype)
df = df[ ,-5]
Loc = rep("", nrow(df))
Loc[df$b < -1] = "Cytosol"
Loc[df$b > 3] = "Nuclear"
df  = df[ ,-1] ## remove L2FC
df = cbind(Loc, df)


set.seed(54321)
inTrain = sample(1:nrow(df), .7*nrow(df) )
train = df[inTrain, ]
test = df[-inTrain, ]

train = train[train$Loc != "", ]
train$Loc = droplevels(train$Loc)

test = test[test$Loc != "", ]
test$Loc = droplevels(test$Loc)


library(caret)
tunegrid <- expand.grid(.mtry=seq(30,48,by = 4))
rf <- train(Loc ~ ., data = train, 
                 method = "rf",
              tuneGrid = tunegrid,
                 ntree = 100)


# mtry  Accuracy   Kappa    
# 30    0.7314787  0.4138936
# 34    0.7325776  0.4167654
# 38    0.7300036  0.4110123
# 42    0.7346429  0.4213275
# 46    0.7325236  0.4160394


test_vals = predict(rf, test, type = "prob")
confusionMatrix( reference = test$Loc, data = test_vals)
saveRDS(rf, "./Data/Raw_RF.rds")

# Confusion Matrix and Statistics
# 
# Reference
# Prediction Cytosol Nuclear
# Cytosol     168      50
# Nuclear     149     438
# 
# Accuracy : 0.7528          
# 95% CI : (0.7215, 0.7822)
# No Information Rate : 0.6062          
# P-Value [Acc > NIR] : < 2.2e-16       
# 
# Kappa : 0.4523          
# Mcnemar's Test P-Value : 3.73e-12        
#                                           
#             Sensitivity : 0.5300          
#             Specificity : 0.8975          
#          Pos Pred Value : 0.7706          
#          Neg Pred Value : 0.7462          
#               Precision : 0.7706          
#                  Recall : 0.5300          
#                      F1 : 0.6280          
#              Prevalence : 0.3938          
#          Detection Rate : 0.2087          
#    Detection Prevalence : 0.2708          
#       Balanced Accuracy : 0.7138 

preds2 = prediction(test_vals[ ,2], test$Loc)
perf2 = performance(preds2, "auc")
