

load("./Data/RFmodel.RData")


trains = read.csv("./Data/Training_preds.csv")
tests = read.csv("./Data/Tests_preds.csv")

trains$Raw_RF = predict(rf, trainRAW)
tests$Raw_RF = predict(rf, testRAW)



lmFit = lm(L2FC ~ DNN + SDA_RF + Raw_RF, data = trains)

vals = predict(lmFit, newdata = tests)

cor.test(vals, tests$L2FC)

plot(tests$Raw_RF, tests$L2FC, xlab = "Predicted", ylab = "Observed", main = "LncRNA Nuclear/Cytosolic ratio")

lm(tests$L2FC ~  tests$Raw_RF)

abline(a = -0.1363, b = 1.1233, lwd = 3, col = "blue")
text(x= -1, y = 2, labels = "Corr =  0.36", cex =1.5)


