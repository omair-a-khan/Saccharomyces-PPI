####################
# Data Preparation #
####################

setwd("~/Dropbox/Documents/School/Spring 2015/STAT 897D")
load("Project2.RData")

Y.train <- Network.train[lower.tri(Network.train)] # 7381 pairs in training # data
n.train <- length(Y.train)
Y.train.mean <- mean(Y.train) # 0.02669015 proportion of PPI interactions in training data

Y.valid <- Network.valid[lower.tri(Network.valid)] # 1378 pairs in validation # data
n.valid <- length(Y.valid)
Y.valid.mean <- mean(Y.valid) # 0.03120464 proportion of PPI interactions in validation data

X.train = NULL # 7381 by 5 matrix with 5 predictors for each gene pair
for (i in 1:(dim(DATA.train)[1]-1))
  for (j in (i+1):dim(DATA.train)[1])
    X.train = rbind(X.train,
                    c(mean(DATA.train[i,]), mean(DATA.train[j,]),
                      cov(DATA.train[i,],DATA.train[j,]),
                      var(DATA.train[i,]), var(DATA.train[j,])))
data.train <- as.data.frame(cbind(Y.train, X.train))
names(data.train) <- c("Y", "X1", "X2", "X3", "X4", "X5")

#  Y: PPI(i,j) **coded as 0 and 1 
# X1: mean of gene i
# X2: mean of gene j
# X3: covariance(i,j)
# X4: variance of gene i
# X5: variance of gene j

X.valid = NULL # 1378 by 5 matrix with 5 predictors for each gene pair
for (i in 1:(dim(DATA.valid)[1]-1))
  for (j in (i+1):dim(DATA.valid)[1])
    X.valid = rbind(X.valid,
              c(mean(DATA.valid[i,]), mean(DATA.valid[j,]),
                cov(DATA.valid[i,],DATA.valid[j,]),
                var(DATA.valid[i,]), var(DATA.valid[j,])))
data.valid <- as.data.frame(X.valid)
names(data.valid) <- c("X1", "X2", "X3", "X4", "X5")

#######################
# Logistic Regression #
#######################

fit.logistic <- glm(Y ~ ., data=data.train, family=binomial("logit"))
summary(fit.logistic)

post.train.logistic <- fit.logistic$fitted.values # n.train posterior probabilities of Y=1
cutoff.logistic <- sort(post.train.logistic, decreasing=T)[201]
#      5272 
# 0.1720898
Ghat.train.logistic <- ifelse(post.train.logistic > cutoff.logistic, 1, 0)
table(Ghat.train.logistic,Y.train) # classification table
#      0    1
# 0 7062  119
# 1  122   78
sum(abs(Ghat.train.logistic-Y.train))/n.train # training data classification error rate = (122+119)/7381
# 0.0326514
sum(Ghat.train.logistic==1&Y.train==1)/sum(Y.train==1) # sensitivity = 78/(119+78)
# 0.3959391
sum(Ghat.train.logistic==0&Y.train==0)/sum(Y.train==0) # specificity = 7062/(7062+122)
# 0.9830178

post.valid.logistic <- predict(fit.logistic, data.valid, type="response") # n.valid post probs
Ghat.valid.logistic <- ifelse(post.valid.logistic > cutoff.logistic,1,0) # use same probability cutoff
table(Ghat.valid.logistic,Y.valid) # classification table
#      0    1
# 0 1301   19
# 1   34   24
sum(abs(Ghat.valid.logistic-Y.valid))/n.valid # classification error rate = (34+19)/1378
# 0.03846154
sum(Ghat.valid.logistic==1&Y.valid==1)/sum(Y.valid==1) # sensitivity = 24/(19+24)
# 0.5581395
sum(Ghat.valid.logistic==0&Y.valid==0)/sum(Y.valid==0) # specificity = 1301/(1301+34)
# 0.9745318

################################
# Linear Discriminant Analysis #
################################

library(MASS)
model.lda <- lda(Y~., data=data.train)
model.lda
plot(model.lda)

post.train.lda <- predict(model.lda)$posterior[,2] # n.train posterior probabilities of Y=1
cutoff.lda <- sort(post.train.lda, decreasing=T)[201]
#      1614 
# 0.1449664 
Ghat.train.lda <- ifelse(post.train.lda>cutoff.lda,1,0) # classification rule
table(Ghat.train.lda,Y.train) # classification table
#      0    1
# 0 7051  130
# 1  133   67
sum(abs(Ghat.train.lda-Y.train))/n.train # training data classification error rate = (133+130)/7381
# 0.03563203
sum(Ghat.train.lda==1&Y.train==1)/sum(Y.train==1) # sensitivity = 67/(130+67)
# 0.3401015
sum(Ghat.train.lda==0&Y.train==0)/sum(Y.train==0) # specificity = 7051/(7051+133)
# 0.9814866

post.valid.lda <- predict(model.lda, data.valid)$posterior[,2] # n.valid posterior probabilities of Y=1
Ghat.valid.lda <- ifelse(post.valid.lda>cutoff.lda,1,0) # use same probability cutoff
table(Ghat.valid.lda,Y.valid) # classification table
#      0    1
# 0 1311   21
# 1   24   22
sum(abs(Ghat.valid.lda-Y.valid))/n.valid # classification error rate = (24+21)/1378
# 0.03265602
sum(Ghat.valid.lda==1&Y.valid==1)/sum(Y.valid==1) # sensitivity = 22/(21+22)
# 0.5116279
sum(Ghat.valid.lda==0&Y.valid==0)/sum(Y.valid==0) # specificity = 1311/(1311+24)
# 0.9820225

###################################
# Quadratic Discriminant Analysis #
###################################

library(MASS)
model.qda <- qda(Y~., data=data.train)
model.qda

post.train.qda <- predict(model.qda)$posterior[,2] # n.train posterior probabilities of Y=1
cutoff.qda <- sort(post.train.qda, decreasing=T)[201]
#      3458 
# 0.3514358 
Ghat.train.qda <- ifelse(post.train.qda>cutoff.qda,1,0) # classification rule
table(Ghat.train.qda,Y.train) # classification table
#      0    1
# 0 7055  126
# 1  129   71
sum(abs(Ghat.train.qda-Y.train))/n.train # training data classification error rate = (129+126)/7381
# 0.03454816
sum(Ghat.train.qda==1&Y.train==1)/sum(Y.train==1) # sensitivity = 71/(126+71)
# 0.3604061
sum(Ghat.train.qda==0&Y.train==0)/sum(Y.train==0) # specificity = 7055/(7055+129)
# 0.9820434

post.valid.qda <- predict(model.qda, data.valid)$posterior[,2] # n.valid posterior probabilities of Y=1
Ghat.valid.qda <- ifelse(post.valid.qda>cutoff.qda,1,0) # use same probability cutoff
table(Ghat.valid.qda,Y.valid) # classification table
#      0    1
# 0 1307   19
# 1   28   24
sum(abs(Ghat.valid.qda-Y.valid))/n.valid # classification error rate = (28+19)/1378
# 0.0341074
sum(Ghat.valid.qda==1&Y.valid==1)/sum(Y.valid==1) # sensitivity = 24/(19+24)
# 0.5581395
sum(Ghat.valid.qda==0&Y.valid==0)/sum(Y.valid==0) # specificity = 1307/(1307+28)
# 0.9790262

#######################
# K Nearest Neighbors #
#######################

X <- rbind(X.train,X.valid)
X.std <- scale(X)
X.train.std <- X.std[1:n.train,]
X.valid.std <- X.std[(n.train+1):(n.train+n.valid),]

library(class)
mer <- rep(NA, 30) # misclassification error rates based on leave-one-out cross-validation

set.seed(2014) # seed must be set because R randomly breaks ties

for (i in 1:30) mer[i] <- sum((Y.train-(c(knn.cv(train=X.train.std, cl=Y.train, k=i))-1))^2)/n.train
plot(mer)
which.min(mer) # minimum occurs at k=13
set.seed(2014)
model.knn <- knn(train=X.train.std, test=X.train.std, cl=Y.train, k=13, prob=T)
predclass.knn <- c(model.knn)-1 # convert factor to numeric classes

predprob.knn <- attr(model.knn, "prob") # proportion of votes for winning class

post.train.knn <- predclass.knn*predprob.knn+(1-predclass.knn)*(1-predprob.knn) # n.train post probs of Y=1

cutoff.knn <- sort(post.train.knn, decreasing=T)[201]
# 0.2307692 

Ghat.train.knn <- ifelse(post.train.knn>cutoff.knn,1,0) # classification rule

table(Ghat.train.knn,Y.train) # classification table
#      0    1
# 0 7091  103
# 1   93   94
sum(abs(Ghat.train.knn-Y.train))/n.train # training data classification error rate = (93+103)/7381
# 0.02655467
sum(Ghat.train.knn==1&Y.train==1)/sum(Y.train==1) # sensitivity = 94/(103+94)
# 0.4771574
sum(Ghat.train.knn==0&Y.train==0)/sum(Y.train==0) # specificity = 7091/(7091+93)
# 0.9870546

set.seed(2014)
model.knn <- knn(train=X.train.std, test=X.valid.std, cl=Y.train, k=13, prob=T)
predclass.knn <- c(model.knn)-1 # convert factor to numeric classes

predprob.knn <- attr(model.knn, "prob") # proportion of votes for winning class

post.valid.knn <- predclass.knn*predprob.knn+(1-predclass.knn)*(1-predprob.knn) # n.valid post probs of Y=1

Ghat.valid.knn <- ifelse(post.valid.knn>cutoff.knn,1,0) # use same probability cutoff

table(Ghat.valid.knn,Y.valid) # classification table
#      0    1
# 0 1294   30
# 1   41   13
sum(abs(Ghat.valid.knn-Y.valid))/n.valid # classification error rate = (41+30)/1378
# 0.05152395
sum(Ghat.valid.knn==1&Y.valid==1)/sum(Y.valid==1) # sensitivity = 13/(30+13)
# 0.3023256
sum(Ghat.valid.knn==0&Y.valid==0)/sum(Y.valid==0) # specificity = 1294/(1294+41)
# 0.9692884

###########################
# Logistic Regression GAM #
###########################

library(gam)

model.gam <- gam(Y ~ s(X1,df=5) + s(X2,df=5) + s(X3,df=5) + s(X4,df=5) + s(X5,df=5), data.train, family=binomial)
summary(model.gam)

par(mfrow=c(1,5))
plot(model.gam, se = TRUE, col="blue")

post.train.gam <- model.gam$fitted.values # n.train posterior probabilities of Y=1

cutoff.gam <- sort(post.train.gam, decreasing=T)[201]
#      3868 
# 0.1725854 
Ghat.train.gam <- ifelse(post.train.gam>cutoff.gam,1,0) # classification rule

table(Ghat.train.gam,Y.train) # classification table
#      0    1
# 0 7071  110
# 1  113   87
sum(abs(Ghat.train.gam-Y.train))/n.train # training data classification error rate = (113+110)/7381
# 0.03021271
sum(Ghat.train.gam==1&Y.train==1)/sum(Y.train==1) # sensitivity = 87/(110+87)
# 0.4416244
sum(Ghat.train.gam==0&Y.train==0)/sum(Y.train==0) # specificity = 7071/(7071+113)
# 0.9842706

post.valid.gam <- predict(model.gam, data.valid, type="response") # n.valid post probs
Ghat.valid.gam <- ifelse(post.valid.gam>cutoff.gam,1,0) # use same probability cutoff
table(Ghat.valid.gam,Y.valid) # classification table
#      0    1
# 0 1254   19
# 1   81   24
sum(abs(Ghat.valid.gam-Y.valid))/n.valid # classification error rate = (81+19)/1378
# 0.07256894
sum(Ghat.valid.gam==1&Y.valid==1)/sum(Y.valid==1) # sensitivity = 19/(19+24)
# 0.5581395
sum(Ghat.valid.gam==0&Y.valid==0)/sum(Y.valid==0) # specificity = 1254/(1254+81)
# 0.9393258

####################################
# Results summary for training set #
####################################

# Method:   Cutoff Probability, Error rate, Sensitivity, Specificity
# Logistic: 17.21%, 3.265%, 39.59%, 98.30%
# LDA:      14.50%, 3.563%, 34.01%, 98.15%
# QDA:      35.14%, 3.455%, 36.04%, 98.20%
# kNN:      23.08%, 2.655%, 47.72%, 98.71%
# GAM:      17.26%, 3.021%, 44.16%, 98.43%

###################################
# Results summary for testing set #
###################################

# Method:   Cutoff Probability, Error rate, Sensitivity, Specificity
# Logistic: 17.21%, 3.845%, 55.81%, 97.45%
# LDA:      14.50%, 3.266%, 51.16%, 98.20%
# QDA:      35.14%, 3.411%, 55.81%, 97.90%
# kNN:      23.08%, 5.152%, 30.23%, 96.93%
# GAM:      17.26%, 7.257%, 55.81%, 93.93%

Error = c(3.845, 3.266, 3.411, 5.152, 7.257)
Sensitivity = c(55.81, 51.16, 55.81, 30.23, 55.81)
Specificity = c(97.45, 98.20, 97.90, 96.93, 93.93)
results <-data.frame(Error, Sensitivity, Specificity)
