# Deep Learning for lncRNA localization
# 2/27/16

library(h2o)
set.seed(54321)
df = read.csv("./Data/4mers.csv", row.names = 1)
# Import data -------------------------------------------------------------
df = read.csv("./Data/lncRNA_featureset.csv", row.names = 1)
filt = -c(1,2,4,5,6,7,8,13)
df = df[ , filt]

Loc = cbind(local = rep("", nrow(df)), df)
Loc$local = as.character(Loc$local)
Loc$local[Loc$log2FoldChange <= 0] = "Cytosol"
Loc$local[Loc$log2FoldChange >= 2] = "Nuclear"
Loc = Loc[ ,-2]
Loc = Loc[Loc$local != "", ]
Loc$local = as.factor(Loc$local)

# Start H2O cluster -------------------------------------------------------
localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = -1, max_mem_size = "14G")

dat_h2o = as.h2o(df)
splits = h2o.splitFrame(dat_h2o, ratios =  c(0.6, 0.2), seed = 54321)
train = splits[[1]]
valid = splits[[2]]
test = splits[[3]]
# dat_h2o = h2o.importFile(localH2O, path = ".../")

hyper_params <- list(
    activation=c("RectifierWithDropout"),
    hidden=list(  c(3800, 3800) ),
    input_dropout_ratio=c( 0.1),
    l1 = c(1e-6, 1e-4, 0),
    l2 = c(1e-6, 1e-4, 0)
)

search_criteria = list(strategy = "RandomDiscrete", max_models = 50, seed=54321, stopping_rounds=5, stopping_tolerance=1e-2)

dl_random_grid <- h2o.grid(
    algorithm="deeplearning",
    grid_id = "dl_grid_random",
    training_frame=train,
    validation_frame=valid, 
    x=2:ncol(train), 
    y=1,
    epochs= 1,
    stopping_metric="RMSE",
    stopping_tolerance=1e-2,        ## stop when RMSE does not improve by >=1% for 5 scoring events
    stopping_rounds=5,            
    hyper_params = hyper_params,
    search_criteria = search_criteria
)          

save.image("./Data/Grid_search.RData")


grid <- h2o.getGrid("dl_grid_random" ,sort_by="RMSE",decreasing=FALSE)
grid
best_model <- h2o.getModel(grid@model_ids[[1]])
# 6/6/17
# best model
# Hyper-Parameter Search Summary: ordered by increasing residual_deviance
# activation               hidden input_dropout_ratio     l1     l2               model_ids  residual_deviance
# 1      TanhWithDropout      [128, 128, 128]                 0.1 5.0E-5 8.0E-5  dl_grid_random_model_6 2.0823344552451615
# 
# 6/8/17
# Model Details:
#     ==============
#     
#     H2ORegressionModel: deeplearning
# Model ID:  dl_grid_random_model_11 
# Status of Neuron Layers: predicting log2FoldChange, regression, gaussian distribution, Quadratic loss, 271,105 weights/biases, 3.8 MB, 23,820 training samples, mini-batch size 1
# layer units             type dropout       l1       l2 mean_rate rate_rms momentum mean_weight weight_rms mean_bias bias_rms
# 1     1  4104            Input 20.00 %                                                                                        
# 2     2    64 RectifierDropout 50.00 % 0.000001 0.004501  0.006798 0.015517 0.000000    0.000120   0.025633  0.358726 0.054511
# 3     3    64 RectifierDropout 50.00 % 0.000001 0.004501  0.001860 0.000796 0.000000   -0.004417   0.122234  0.886999 0.033358
# 4     4    64 RectifierDropout 50.00 % 0.000001 0.004501  0.004299 0.005892 0.000000   -0.019839   0.117244  0.829142 0.133902
# 5     5     1           Linear         0.000001 0.004501  0.000231 0.000074 0.000000   -0.011297   0.093161  0.054882 0.000000





test_vals = h2o.predict(object = best_model, test)
plot(as.vector(as.numeric(test_vals)), as.vector(as.numeric(test$log2FoldChange) ))
abline(lm(as.matrix(test_vals) ~ as.matrix(test$log2FoldChange)), col="blue", lwd=2)
cor.test(as.vector(as.numeric(test_vals)), as.vector(as.numeric(test$log2FoldChange) ))

H2O = h2o.deeplearning(x =2:ncol(train), #predictors
                       y = 1, #label 
                       training_frame = train,
                       activation = "RectifierWithDropout",
                       input_dropout_ratio = 0.2,
                       validation_frame = valid,
                       stopping_rounds = 10,
                       stopping_metric = "RMSE",
                       stopping_tolerance = 0.01,
                       seed = 54321,
                       hidden = c(500,500,500),
                       epochs = 1000)


test_vals = h2o.predict(object = H2O, test)
plot(as.vector(as.numeric(test_vals)), as.vector(as.numeric(test$log2FoldChange) ), ylim = range(-1,4), xlim = range(-1, 4))
abline(lm(as.matrix(test_vals) ~ as.matrix(test$log2FoldChange)), col="blue", lwd=2)
cor.test(as.vector(as.numeric(test_vals)), as.vector(as.numeric(test$log2FoldChange) ))


hyper_params <- list(
    activation=c("RectifierWithDropout","Tanh","Rectifier"),
    hidden=list(c(128,128, 128),c(500,500,500),c(2000,2000,2000), c(128,128,128,128),c(500,500,500,500)),
    input_dropout_ratio=c(0,0.1,0.2),
    l1=seq(0,1e-4,1e-6),
    l2=seq(0,1e-4,1e-6)
)

search_criteria = list(strategy = "RandomDiscrete", max_models = 100, seed=54321, stopping_rounds=5, stopping_tolerance=1e-2)

dl_random_grid <- h2o.grid(
    algorithm="deeplearning",
    grid_id = "dl_grid_random",
    training_frame=train,
    validation_frame=as.h2o(validation), 
    x=2:4101, 
    y=1,
    epochs= 10,
    stopping_metric="RMSE",
    stopping_tolerance=1e-2,        ## stop when RMSE does not improve by >=1% for 2 scoring events
    stopping_rounds=5,  
    mini_batch_size = 5,
    hyper_params = hyper_params,
    search_criteria = search_criteria
)                                



grid <- h2o.getGrid("dl_grid_random",sort_by="RMSE",decreasing=FALSE)
grid

best_model <- h2o.getModel(grid@model_ids[[1]])
valid_vals = h2o.predict(object = best_model, as.h2o(validation))
cor(as.vector(valid_vals), validation$log2FoldChange)
# 
# 
# Model Details:
#     ==============
#     
#     H2ORegressionModel: deeplearning
# Model ID:  dl_grid_random_model_21 
# Status of Neuron Layers: predicting log2FoldChange, regression, gaussian distribution, Quadratic loss, 575,233 weights/biases, 7.3 MB, 43,325 training samples, mini-batch size 5
# layer units   type dropout       l1       l2 mean_rate rate_rms momentum mean_weight weight_rms mean_bias bias_rms
# 1     1  4105  Input 20.00 %                                                                                        
# 2     2   128   Tanh  0.00 % 0.000000 0.000088  0.136998 0.123641 0.000000   -0.000174   0.030224  0.007384 0.060973
# 3     3   128   Tanh  0.00 % 0.000000 0.000088  0.081193 0.037227 0.000000   -0.000105   0.080598 -0.006018 0.067926
# 4     4   128   Tanh  0.00 % 0.000000 0.000088  0.129371 0.056110 0.000000    0.000183   0.068002  0.007204 0.078359
# 5     5   128   Tanh  0.00 % 0.000000 0.000088  0.402972 0.386623 0.000000    0.000221   0.060661  0.002537 0.061877
# 6     6     1 Linear         0.000000 0.000088  0.009211 0.009362 0.000000   -0.006868   0.076838  0.014995 0.000000
# RMSE:  1.17127



hyper_params <- list(
    activation=c("Tanh","TanhWithDropout"),
    hidden=list( c(128,128,128,128),c(500,500,500,500), c(1000,1000,1000,1000), c(128,128,128,128,128)),
    input_dropout_ratio=c(0,0.1,0.2)
)

search_criteria = list(strategy = "RandomDiscrete", max_models = 100, seed=54321, stopping_rounds=5, stopping_tolerance=1e-2)

dl_random_grid <- h2o.grid(
    algorithm="deeplearning",
    grid_id = "dl_grid_random",
    training_frame=train,
    validation_frame=as.h2o(validation), 
    x=2:4101, 
    y=1,
    epochs= 10,
    l2 = 0.000088,
    stopping_metric="RMSE",
    stopping_tolerance=1e-2,        ## stop when RMSE does not improve by >=1% for 2 scoring events
    stopping_rounds=5,  
    mini_batch_size = 1,
    hyper_params = hyper_params,
    search_criteria = search_criteria
)                                

grid <- h2o.getGrid("dl_grid_random",sort_by="RMSE",decreasing=FALSE)
grid

best_model <- h2o.getModel(grid@model_ids[[1]])
best_model

# Model Details:
#     ==============
#     
#     H2ORegressionModel: deeplearning
# Model ID:  dl_grid_random_model_11 
# Status of Neuron Layers: predicting log2FoldChange, regression, gaussian distribution, Quadratic loss, 575,233 weights/biases, 7.3 MB, 62,376 training samples, mini-batch size 1
# layer units   type dropout       l1       l2 mean_rate rate_rms momentum mean_weight weight_rms mean_bias bias_rms
# 1     1  4105  Input 10.00 %                                                                                        
# 2     2   128   Tanh  0.00 % 0.000000 0.000088  0.040847 0.045541 0.000000    0.000405   0.038455 -0.008847 0.055903
# 3     3   128   Tanh  0.00 % 0.000000 0.000088  0.021916 0.011745 0.000000    0.001181   0.101030  0.006860 0.066021
# 4     4   128   Tanh  0.00 % 0.000000 0.000088  0.037036 0.013879 0.000000    0.000252   0.092952 -0.005113 0.090580
# 5     5   128   Tanh  0.00 % 0.000000 0.000088  0.205046 0.191894 0.000000    0.000044   0.068054  0.000786 0.078641
# 6     6     1 Linear         0.000000 0.000088  0.001894 0.001398 0.000000   -0.005971   0.065559 -0.063718 0.000000
# 
# 
# H2ORegressionMetrics: deeplearning
# ** Reported on training data. **
#     ** Metrics reported on full training frame **
#     
#     MSE:  0.2586711
# RMSE:  0.5085972
# MAE:  0.3647542
# RMSLE:  NaN
# Mean Residual Deviance :  0.2586711
# 
# 
# H2ORegressionMetrics: deeplearning
# ** Reported on validation data. **
#     ** Metrics reported on full validation frame **
#     
#     MSE:  1.150288
# RMSE:  1.072515
# MAE:  0.7103557
# RMSLE:  NaN
# Mean Residual Deviance :  1.150288

valid_vals = h2o.predict(object = best_model, as.h2o(validation))
cor(as.vector(valid_vals), validation$log2FoldChange)^2
sqrt(mean((as.vector(valid_vals) - validation$log2FoldChange)^2 ))
plot(as.vector(valid_vals), validation$log2FoldChange, xlab = "Prediction",ylab = "Observed", main = "Validation Set: LncRNA Subcellular Log2 Fold-change")
legend(x = "bottomright", legend = paste("Pearson Correlation =", round(cor(as.vector(valid_vals), validation$log2FoldChange),3)))




H2O = h2o.deeplearning(x =2:4101, #predictors
                       y = 1, #label 
                       training_frame = dat_h2o,
                       activation = "Tanh",
                       input_dropout_ratio = 0.1,
                       stopping_rounds = 15,
                       stopping_metric = "RMSE",
                       stopping_tolerance = 0.0005,
                       seed = 54321,
                       l1 = 0,
                       l2 = 0.000088,
                       hidden = c(128,128,128,128),
                       epochs = 200,
                       variable_importances = TRUE)
save.image("DeepLearn_Final.RData")

test_vals = h2o.predict(object = H2O, test)
cor(as.vector(test_vals), as.vector(test$log2FoldChange))

pdf("./Figures/Test_set_error.pdf")
par(mar=c(4.5,5,3,1))
plot(x= as.vector(test$log2FoldChange), y = as.vector(test_vals), xlab = "Observed Log2 Fold-Change", ylab = "Expected Log2 Fold-Change", main = "Test set: Observed vs Predicted", cex.lab=2.5, cex.main =2)
abline(lm(as.matrix(valid_vals) ~ as.matrix(test$log2FoldChange)), col="blue", lwd=2)
text(x=0, y= -0.8, labels = paste("Pearson Correlation =", round(cor(as.vector(test_vals), as.vector(test$log2FoldChange)),3)), cex=2)
text(x=0, y= -1.2, labels = paste("RMSE = ",RMSE), cex=2)
dev.off()


RMSE=round(sqrt(mean((as.vector(test_vals) - as.vector(test$log2FoldChange))^2 )),3)
var_imp = h2o.varimp(H2O)


pdf("./Figures/Variable_importance.pdf")
par(mar=c(9,4,2,1))
barplot(var_imp[1:20,4]*100, names.arg = var_imp[1:20,1], las=2, col ="blue2", main = "Relative Variable Importance", cex.axis = 1.5, cex.names = 1.5, cex.main=2)
dev.off()

write.csv(var_imp, "./Data/Variable_imp.csv")

for ( i in 4:50){
    var = var_imp$variable[i]
    #print(var_imp$variable[i])
   # print(cor(var_imp$))
    print(cor(df[colnames(df) == var], df$log2FoldChange))
}
    
