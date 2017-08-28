
# RF model ----------------------------------------------------------------
library(h2o)
df = read.csv("./Data/4mers.csv", row.names = 1)

localH2O = h2o.init(ip="localhost", port = 54321, startH2O = TRUE, nthreads = 2, max_mem_size = "10G")

dat_h2o = as.h2o(df)
# train_un = dat_h2o[Loc_un, ]
# 
# dat_h2o = dat_h2o[-Loc_un, ]
# 
# dat_h2o = h2o.cbind(dat_h2o, as.h2o(droplevels(resp[ - Loc_un])))

splits = h2o.splitFrame(dat_h2o, ratios =  c(0.7), seed = 54321)
train_hex = splits[[1]]
test_hex = splits[[2]]


trainRAW = as.data.frame(train_hex)
testRAW = as.data.frame(test_hex)
h2o.shutdown(FALSE)
rm(dat_h2o, localH2O)


library(randomForest)
rf = randomForest(x = trainRAW[ ,-1],
                  y = trainRAW$log2FoldChange,
                  data = trainRAW,
                  ntree=500,
                  do.trace = TRUE)

levels(testRAW$biotype) = levels(trainRAW$biotype)

test_vals = predict(rf, newdata = testRAW)

test_vals = predict(rf, x = testRAW[,-1], y = testRAW$log2FoldChange)

sqrt(mean((test_vals - testRAW$log2FoldChange)^2))
plot(test_vals, testRAW$log2FoldChange)
abline(lm(testRAW$log2FoldChange ~ test_vals), lwd =2, col= "blue")

save.image("./Data/RFmodel.RData")
