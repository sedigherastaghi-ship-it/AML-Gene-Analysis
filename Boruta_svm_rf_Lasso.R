library(Hmisc)
library(StepReg)
library(pROC)
library(ggplot2)
###########Feature selection for machine learning##########
library(dplyr)
library(caret)
#nzv <- nearZeroVar(x3)
library(Boruta)

########3Feature selection####

##1- data prparation
#data <- read.delim("E:/papers/MLS/luk/Data& Codes/GSE9476_series_matrix-data.txt", header=FALSE)
data <- read.delim(file.choose(), header=FALSE)
x=as.matrix(t(data[,-1]))
sum(is.na(x))
colnames(x)=as.vector(data[,1])
y=factor(c(rep(0,18),rep(1,26), rep(0,20))) 
data.lukemia=data.frame(y,x)

## subset selcetion in data
y3=y[19:64]
x3=x[c(19:64),]
y3=factor(y3, labels = c("Normal","AML"))
table(y3)

data3=data.frame(x3,y3)
######Boruta Feature Elimination
set.seed(2026)
boruta_res <- Boruta(y3 ~ .,pValue = 0.0005, data = data3, doTrace = 2, ntree = 500)
print(boruta_res)
par(mfrow=c(1,1))
plot(boruta_res, xlab = "", xaxt = "n", main = "Boruta Feature Selection")
#final_vars <- getSelectedAttributes(boruta_res, withTentative = FALSE)
final_vars=c("X200736_s_at" ,"X203195_s_at" ,"X203556_at" ,  "X204949_at" ,  "X205801_s_at", "X207992_s_at", "X208772_at"  ,
             "X208773_s_at", "X209112_at" ,  "X209586_s_at", "X210183_x_at" ,"X212184_s_at", "X212549_at" ,  "X212560_at",  
             "X213351_s_at" ,"X213506_at" ,  "X213511_s_at" ,"X215785_s_at", "X216095_x_at", "X217756_x_at", "X218999_at" , 
             "X219315_s_at", "X219966_x_at", "X221556_at" ,  "X221771_s_at", "X221979_at" ,  "X222307_at"  , "X34858_at" ,  
             "X78383_at")
data3_new <- data3[, c(final_vars, "y3")]
xk.boruta=as.matrix(data3[, c(final_vars)])
a=colnames(xk.boruta)

######LASSO#############
library(glmnet)
cvfit=cv.glmnet(xk.boruta,y3, family="binomial", alpha=1,nfolds=10 )
preds <- predict(cvfit, newx=as.matrix(xk.boruta), s="lambda.min", type="response")
library(pROC)
#roc.lasso=roc(y3,as.vector(preds),plot=TRUE,print.auc=TRUE,main='TOTAL',col="blue3",lty=1,size=3)
auc3=auc(y3,as.vector(preds))
b.total=as.matrix(coef(cvfit,s="lambda.min")[-1,]) 
b.sig=(b.total!=0) ;sig.b=b.total[b.sig==1,]
cbind(names(sig.b),as.vector(sig.b),exp(as.vector(sig.b)) )

t.sys0=Sys.time()
y.sim=y3
t1=table(y.sim)[2];t0=table(y.sim)[1]
nsim=1000
auc.lasso=auc.ens=0
b.lasso=matrix(0,nrow=ncol(xk.boruta)+1,ncol=nsim)

for(bs in 1:nsim){
  id.sample1=sample(1:t1,t1,replace =T)
  id.sample0=sample((t1+1):(t1+t0),t0,replace =T)
  id.sample=c(id.sample1,id.sample0)
  y3=y.sim[id.sample]
  xkmem=xk.boruta[id.sample,]
  ## logitic regression by elasticnet approach
  cvfit=cv.glmnet(xkmem,y3, family="binomial", alpha=1,nfolds=10 )
  lambda_best <- cvfit$lambda.min
  preds <- predict(cvfit, newx=as.matrix(xkmem), s=lambda_best, type="response")
  #pROC::roc(y3,as.vector(preds1),plot=TRUE,print.auc=TRUE,main='CLUSTER=1',col="blue3")
  auc.lasso[bs]=auc(y3,as.vector(preds))
  
  b.lasso[,bs]=coef(cvfit,s=lambda_best)[,1]

}

b.LASSO=t(b.lasso[-1,])
colnames(b.LASSO)=a

b.LASSOi=se.b.LASSOi=lcibeta=ucibeta=oddsi=lciodds=uciodds=0
for(i in 1:ncol(b.LASSO)){
  b.LASSOi[i]=mean(b.LASSO[,i])
  se.b.LASSOi[i]=sd(b.LASSO[,i])/sqrt(nrow(x3))
  lcibeta[i]=b.LASSOi[i]+qnorm(0.00025)*se.b.LASSOi[i]
  ucibeta[i]=b.LASSOi[i]+qnorm(0.99975)*se.b.LASSOi[i]
  oddsi[i]=exp(b.LASSOi[i])
  lciodds[i]=exp(lcibeta[i])
  uciodds[i]=exp(ucibeta[i])
}
table3.lasso=round(cbind(b.LASSOi,se.b.LASSOi,lcibeta,ucibeta,oddsi,lciodds,uciodds),3)

df_genes.LASSO <- data.frame(
  gene = a,
  mean = b.LASSOi,
  lower = lcibeta,
  upper = ucibeta
)

df_filtered.LASSO <- df_genes.LASSO %>%
  arrange(mean) %>%
  mutate(gene = factor(gene, levels = gene))
IM.LASSO=ggplot(df_filtered.LASSO, aes(x = mean, y = gene)) +
  
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3, size = 0.8,color="#1f77b4") +
  geom_point(size = 2, color = "#1f77b4") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "LASSO",
    x = "Mean of Beta with CI",
    y = "Gene ID"
  ) +
  theme(
    axis.text.y = element_text(size = 8), 
    panel.grid.minor = element_blank())+
  theme_bw() + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black")
  )

############elastic net#############
y3=y[19:64]
x3=x[c(19:64),]
y3=factor(y3, labels = c("Normal","AML"))
x3=data.frame(x3)
xk.boruta=as.matrix(x3[, c(final_vars)])
library(glmnet)
cvfit=cv.glmnet(xk.boruta,y3, family="binomial", alpha=0.5,nfolds=10 )
preds <- predict(cvfit, newx=as.matrix(xk.boruta), s="lambda.min", type="response")
library(pROC)
roc.elnet=roc(y3,as.vector(preds),plot=TRUE,print.auc=TRUE,main='TOTAL',col="blue3",,lty=1,size=3)
auc3=auc(y3,as.vector(preds))
b.total=as.matrix(coef(cvfit,s="lambda.min")[-1,]) 
b.sig=(b.total!=0) ;sig.b=b.total[b.sig==1,]
cbind(names(sig.b),as.vector(sig.b),exp(as.vector(sig.b)) )

y.sim=y3
t1=table(y.sim)[2];t0=table(y.sim)[1]
nsim=1000
auc.elnet=auc.ens=0
b.elnet=matrix(0,nrow=ncol(xk.boruta)+1,ncol=nsim)

for(bs in 1:nsim){
  id.sample1=sample(1:t1,t1,replace =T)
  id.sample0=sample((t1+1):(t1+t0),t0,replace =T)
  id.sample=c(id.sample1,id.sample0)
  y3=y.sim[id.sample]
  xkmem=xk.boruta[id.sample,]
  ## logitic regression by elasticnet approach
  cvfit=cv.glmnet(xkmem,y3, family="binomial", alpha=0.5,nfolds=10 )
  lambda_best <- cvfit$lambda.min
  preds <- predict(cvfit, newx=as.matrix(xkmem), s=lambda_best, type="response")
  #pROC::roc(y3,as.vector(preds1),plot=TRUE,print.auc=TRUE,main='CLUSTER=1',col="blue3")
  auc.elnet[bs]=auc(y3,as.vector(preds))
  b.elnet[,bs]=coef(cvfit,s=lambda_best)[,1]
}
a=colnames(xk.boruta)
b.elnet=t(b.elnet[-1,])
colnames(b.elnet)=a

b.elneti=se.b.elneti=lcibeta=ucibeta=oddsi=lciodds=uciodds=0
for(i in 1:ncol(b.elnet)){
  b.elneti[i]=mean(b.elnet[,i])
  se.b.elneti[i]=sd(b.elnet[,i])/sqrt(nrow(x3))
  lcibeta[i]=b.elneti[i]+qnorm(0.00025)*se.b.elneti[i]
  ucibeta[i]=b.elneti[i]+qnorm(0.99975)*se.b.elneti[i]
  oddsi[i]=exp(b.elneti[i])
  lciodds[i]=exp(lcibeta[i])
  uciodds[i]=exp(ucibeta[i])
}
table3.elnet=round(cbind(b.elneti,se.b.elneti,lcibeta,ucibeta,oddsi,lciodds,uciodds),3)

df_genes <- data.frame(
  gene = a,
  mean = b.elneti,
  lower = lcibeta,
  upper = ucibeta
)

df_filtered <- df_genes %>%
  arrange(mean) %>%
  mutate(gene = factor(gene, levels = gene))

IM.en=ggplot(df_filtered, aes(x = mean, y = gene)) +
  
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3, size = 0.8,color="#1f77b4") +
  geom_point(size = 2, color = "#1f77b4") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Elastic Net",
    x = "Mean of Beta with CI",
    y = "Gene ID"
  ) +
  theme(
    axis.text.y = element_text(size = 8), 
    panel.grid.minor = element_blank())+
  theme_bw() + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black")
  )
#############Random Forest######
library(randomForest)
y3=y[19:64]
x3=x[c(19:64),]
y3=factor(y3, labels = c("Normal","AML"))
x3=data.frame(x3)
xk.boruta=as.matrix(x3[, c(final_vars)])
y_final <- as.factor(y3)
mtry_value <- floor(sqrt(ncol(xk.boruta)))

rf_model <- randomForest(x = xk.boruta, 
                         y = factor(y3), 
                         ntree = 1000, 
                         importance = TRUE, 
                         mtry = mtry_value)
#print(rf_model)
rf_probs <- rf_model$votes
prob_class_rf <- attr(rf_probs, "probabilities")[, 2]
importance.rf=(rf_model$importance[,3]/max(rf_model$importance[,3]))*100

y.sim=y3
t1=table(y.sim)[2];t0=table(y.sim)[1]
nsim=1000
auc.rf=auc.ens=0
importance.rf=matrix(0,nrow=ncol(xk.boruta),ncol=nsim)

for(bs in 1:nsim){
  id.sample1=sample(1:t1,t1,replace =T)
  id.sample0=sample((t1+1):(t1+t0),t0,replace =T)
  id.sample=c(id.sample1,id.sample0)
  y3=y.sim[id.sample]
  xk.boruta=xk.boruta[id.sample,]
  rf_model <- randomForest(x = xk.boruta, 
                           y = factor(y3), 
                           ntree = 1000, 
                           importance = TRUE, 
                           mtry = mtry_value)
 
  rf_probs <- rf_model$votes[, 2]
  
  auc.rf[bs]=auc(y3,as.vector(rf_probs))
  importance.rf[,bs]=(rf_model$importance[,3]/max(rf_model$importance[,3]))*100
  
}
a=colnames(xk.boruta)
importance.rf=t(importance.rf[,])
colnames(importance.rf)=a

importance.rfi=se.importance.rfi=lcimport.rf=ucimport.rf=0
for(i in 1:ncol(importance.rf)){
  importance.rfi[i]=mean(importance.rf[,i])
  se.importance.rfi[i]=sd(importance.rf[,i])/sqrt(nrow(x3))
  lcimport.rf[i]=importance.rfi[i]+qnorm(0.00025)*se.importance.rfi[i]
  ucimport.rf[i]=importance.rfi[i]+qnorm(0.99975)*se.importance.rfi[i]
 
  
}
table3.rf=round(cbind(importance.rfi,se.importance.rfi,lcimport.rf,ucimport.rf),3)

rf_genes <- data.frame(
  gene = a,
  mean = importance.rfi,
  lower = lcimport.rf,
  upper = ucimport.rf
)
library(dplyr)
rf_filtered <- rf_genes %>%
  arrange(mean) %>%
  mutate(gene = factor(gene, levels = gene))

IM.rf=ggplot(rf_filtered, aes(x = mean, y = gene)) +
  
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3, size = 0.8,color="#1f77b4") +
  geom_point(size = 2, color = "#1f77b4") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Random Forest",
    x = "Mean of RF Importance with CI",
    y = "Gene ID"
  ) +
  theme(
    axis.text.y = element_text(size = 8), 
    panel.grid.minor = element_blank())+
  theme_bw() + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black")
  )
###########SVM###########
library(e1071)
library(caret)
y3=y[19:64]
x3=x[c(19:64),]
y3=factor(y3, labels = c("Normal","AML"))
x3=data.frame(x3)
xk.boruta=as.matrix(x3[, c(final_vars)])
X_svm <- data.frame(xk.boruta)
y_svm <- as.factor(y3)
train_control <- trainControl(method = "cv", number = 5)
svm_tuned <- train(x = X_svm, 
                   y = y_svm, 
                   method = "svmRadial", 
                   trControl = train_control,
                   preProcess = c("center", "scale"), 
                   tuneLength = 10) 
print("Best SVM Parameters:")
print(svm_tuned$bestTune)
print(svm_tuned)
svm_prob_model <- svm(x = X_svm, 
                      y = y_svm, 
                      kernel = "radial", 
                      probability = TRUE, 
                      scale = TRUE)

svm_probs <- predict(svm_prob_model, X_svm, probability = TRUE)
prob_class_1 <- attr(svm_probs, "probabilities")[, 1]

importance <- varImp(svm_tuned, scale = T) 
print(importance)

y.sim=y3
t1=table(y.sim)[2];t0=table(y.sim)[1]
nsim=1000
auc.svm=auc.ens=0
importance.svm=matrix(0,nrow=ncol(xk.boruta),ncol=nsim)

for(bs in 1:nsim){
  id.sample1=sample(1:t1,t1,replace =T)
  id.sample0=sample((t1+1):(t1+t0),t0,replace =T)
  id.sample=c(id.sample1,id.sample0)
  y.svm=y.sim[id.sample]
  x.svm=data.frame(xk.boruta[id.sample,])
  train_control <- trainControl(method = "cv", number = 5)
  svm_t <- train(x = x.svm, 
                     y = y.svm, 
                     method = "svmRadial", 
                     trControl = train_control,
                     preProcess = c("center", "scale"), 
                     tuneLength = 10)
  importance=varImp(svm_t, scale = T)
  importance.svm[,bs]= importance$importance[,1]
  svm_prob_model <- svm(x = x.svm, 
                        y = y.svm, 
                        kernel = "radial", 
                        probability = TRUE, 
                        scale = TRUE)
  
  svm_probs <- predict(svm_prob_model, x.svm, probability = TRUE)
  svm_probs <- attr(svm_probs, "probabilities")[, 1]
  auc.svm[bs]=auc(y3,as.vector(svm_probs))
}
a=colnames(xk.boruta)
importance.svm=t(importance.svm[,])
colnames(importance.svm)=a

importance.svmi=se.importance.svmi=lcimport.svm=ucimport.svm=0
for(i in 1:ncol(importance.svm)){
  importance.svmi[i]=mean(importance.svm[,i])
  se.importance.svmi[i]=sd(importance.svm[,i])/sqrt(nrow(x3))
  lcimport.svm[i]=importance.svmi[i]+qnorm(0.00025)*se.importance.svmi[i]
  ucimport.svm[i]=importance.svmi[i]+qnorm(0.99975)*se.importance.svmi[i]
 }
table3.svm=round(cbind(importance.svmi,se.importance.svmi,lcimport.svm,ucimport.svm),3)
svm_genes <- data.frame(
  gene = a,
  mean = importance.svmi,
  lower = lcimport.svm,
  upper = ucimport.svm
)
library(dplyr)
svm_filtered <- svm_genes %>%
  arrange(mean) %>%
  mutate(gene = factor(gene, levels = gene))
IM.svm=ggplot(svm_filtered, aes(x = mean, y = gene)) +
   geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.3, size = 0.8,color="#1f77b4") +
  geom_point(size = 2, color = "#1f77b4") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Support Vector Machine",
    subtitle = "",
    x = "Mean of SVM Importance with CI",
    y = "Gene ID"
  ) +
  theme(
    axis.text.y = element_text(size = 8), 
    panel.grid.minor = element_blank())+
  theme_bw() + 
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_line(color = "gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black")
   )
library(patchwork)
print(IM.rf+IM.svm+IM.LASSO+IM.en)
combined<- IM.rf+IM.svm+IM.LASSO+IM.en
print(combined)
ggsave("E://papers/MLS/PLOS one/Figs/Fig8.png", plot = combined, width = 8, height = 10, dpi = 600)
