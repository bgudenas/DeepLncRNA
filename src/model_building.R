# Model building ----------------------------------------------------------
library(caret)
library(randomForest)

df = read.csv("./Data/DE_lncRNAs_kmertable.csv", row.names = 1)

set.seed(54321)
filt = -c(1,3,4,5,6,7)
df = df[ , filt]
inTrain =  sample(1:nrow(df), round(.60*nrow(df)))
training = df[inTrain,  ]
testing = df[ -inTrain,  ]
inValid = sample(1:nrow(testing), round(0.5*nrow(testing)))
validation = testing[inValid, ]
testing = testing[-inValid, ]



RF = train(log2FoldChange ~ ., data = training, method = "rf", ntree = 501)

sqrt(mean((pred - training$log2FoldChange)^2 ))

valid_vals = predict(object = RF, validation)
sqrt(mean((valid_vals - validation$log2FoldChange)^2 ))



rf_caret
2    6.046111  0.4597148
vars = varImp(rf_caret)
rf variable importance



DNN = train(log2FoldChange ~., data= training, method="dnn")
# 
# 
# biotype_mat = model.matrix( ~ biotype - 1, data = training)
# training = training[ , -2]
# training = cbind(biotype_mat ,  training)
# dnn <- sae.dnn.train(as.matrix(training[ ,-6]), training$log2FoldChange , hidden = c(5, 5, 5, 5), sae_output = "sigm", batchsize = 400, numepochs = 5)
# 
# valid_vals = predict(DNN, validation)
# sqrt(mean((valid_vals - validation$log2FoldChange)^2 ))
