# Notes:

library(grf)
if(packageVersion("grf") < '0.10.2') {
  warning("This script requires grf 0.10.2 or higher")
}
library(sandwich)
library(lmtest)
library(Hmisc)
library(ggplot2)

setwd("/Users/BrianChung/Dropbox/Github")

library(haven)
DataQJE <- read.csv("WorkFromHomeProject/treatment_effect.csv")

DataQJE <- subset(DataQJE, year_week!=201049)

DF = DataQJE[,]
person.id = as.numeric(DataQJE$personid)

DF[DF=="yes"]<-1
DF[DF=="no"]<-0
DF$children <- as.numeric(as.character(DF$children))

W = DF$experiment_treatment
Y = DF$perform1
X.raw = cbind(DF[,1:7],DF[,9:9],DF[,11:14],DF[,38:38])
colnames(X.raw)[8] <- "children"
colnames(X.raw)[13] <- "perform10"

year_week.exp = model.matrix(~ factor(X.raw$year_week) + 0)
personid.exp = model.matrix(~ factor(X.raw$personid) + 0)

X = cbind(X.raw[,-which(names(X.raw) %in% c("year_week", "personid"))], year_week.exp, personid.exp)

model1 <-miceadds::lm.cluster(DF,perform1 ~ experiment_treatment+perform10+factor(year_week)+factor(personid),cluster="personid")
summary(model1)


#
# Grow a forest. Add extra trees for the causal forest.
#
  
Y.forest = regression_forest(X, Y, clusters = person.id, equalize.cluster.weights = TRUE)

Y.hat = predict(Y.forest)$predictions
W.forest = regression_forest(X, W, clusters = person.id, equalize.cluster.weights = TRUE)
W.hat = predict(W.forest)$predictions

cf.raw = causal_forest(X, Y, W,
                       Y.hat = Y.hat, W.hat = W.hat,
                       clusters = person.id,
                       equalize.cluster.weights = TRUE)
varimp = variable_importance(cf.raw)
selected.idx = which(varimp > mean(varimp))

cf = causal_forest(X[,selected.idx], Y, W,
                   Y.hat = Y.hat, W.hat = W.hat,
                   clusters = person.id,
                   equalize.cluster.weights = TRUE,
                   tune.parameters = "all")
tau.hat = predict(cf)$predictions

#
# Estimate ATE
#

ATE = average_treatment_effect(cf)
paste("95% CI for the ATE:", round(ATE[1], 3),
      "+/-", round(qnorm(0.975) * ATE[2], 3))

#
# Omnibus tests for heterogeneity
#

# Run best linear predictor analysis
test_calibration(cf)

# Compare regions with high and low estimated CATEs
high_effect = tau.hat > median(tau.hat)
ate.high = average_treatment_effect(cf, subset = high_effect)
ate.low = average_treatment_effect(cf, subset = !high_effect)
paste("95% CI for difference in ATE:",
      round(ate.high[1] - ate.low[1], 3), "+/-",
      round(qnorm(0.975) * sqrt(ate.high[2]^2 + ate.low[2]^2), 3))

#
# Analaysis without fitting the propensity score
#

cf.noprop = causal_forest(X[,selected.idx], Y, W,
                          Y.hat = Y.hat, W.hat = mean(W),
                          tune.parameters = "all",
                          equalize.cluster.weights = TRUE,
                          clusters = person.id)
tau.hat.noprop = predict(cf.noprop)$predictions

ATE.noprop = average_treatment_effect(cf.noprop)
paste("95% CI for the ATE:", round(ATE.noprop[1], 3),
      "+/-", round(qnorm(0.975) * ATE.noprop[2], 3))

pdf("WorkFromHomeProject/tauhat_noprop.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
plot(tau.hat, tau.hat.noprop,
     xlim = range(tau.hat, tau.hat.noprop),
     ylim = range(tau.hat, tau.hat.noprop),
     xlab = "orthogonalized causal forest estimates",
     ylab = "non-orthogonalized causal forest")
abline(0, 1, lwd = 2, lty = 2, col = 4)
par = pardef
dev.off()

#
# Make some graphs
#

pdf("WorkFromHomeProject/tauhat_hist.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
hist(tau.hat, xlab = "estimated CATE", main = "")
dev.off()

pdf("WorkFromHomeProject/tauhat_hist_noprop.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
hist(tau.hat.noprop, xlab = "estimated CATE", main = "")
dev.off()

pdf("WorkFromHomeProject/tauhat_vs_PerformPrev.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
boxplot(tau.hat ~ round(X$perform10), xlab = "Preform Prev", ylab = "estimated CATE")
#lines(smooth.spline(4 + X[,"perform10"], tau.hat, df = 4), lwd = 2, col = 4)
dev.off()

pdf("WorkFromHomeProject/tauhat_vs_Children.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
boxplot(tau.hat ~ round(X$children), xlab = "Have Children", ylab = "estimated CATE")
#lines(smooth.spline(4 + X[,"children"], tau.hat, df = 4), lwd = 2, col = 4)
dev.off()

pdf("WorkFromHomeProject/tauhat_vs_Women.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
boxplot(tau.hat ~ round(X$women), xlab = "Women", ylab = "estimated CATE")
#lines(smooth.spline(4 + X[,"women"], tau.hat, df = 4), lwd = 2, col = 4)
dev.off()

pdf("WorkFromHomeProject/tauhat_vs_LiveWithSpouse.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
boxplot(tau.hat ~ round(X$livewspouse), xlab = "Live With Spouse", ylab = "estimated CATE")
#lines(smooth.spline(4 + X[,"livewspouse"], tau.hat, df = 4), lwd = 2, col = 4)
dev.off()

pdf("WorkFromHomeProject/tauhat_vs_ShortTenure.pdf")
pardef = par(mar = c(5, 4, 4, 2) + 0.5, cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
boxplot(tau.hat ~ round(X$short_tenure), xlab = "Short Tenure", ylab = "estimated CATE")
#lines(smooth.spline(4 + X[,"short_tenure"], tau.hat, df = 4), lwd = 2, col = 4)
dev.off()

#
# Analysis ignoring clusters
#

cf.noclust = causal_forest(X[,selected.idx], Y, W,
                           Y.hat = Y.hat, W.hat = W.hat,
                           tune.parameters = "all")

ATE.noclust = average_treatment_effect(cf.noclust)
paste("95% CI for the ATE:", round(ATE.noclust[1], 3),
      "+/-", round(qnorm(0.975) * ATE.noclust[2], 3))

