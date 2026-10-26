# Version info: R 4.2.2, Biobase 2.58.0, GEOquery 2.66.0, limma 3.54.0
################################################################
#   Data plots for selected GEO samples
library(GEOquery)
library(bigmemory)
gse<- getGEO(filename = "F:\\papers\\MLS\\BMC BIoinformatics\\GSE13159_series_matrix.txt\\GSE13159_series_matrix.txt", getGPL = F, GSEMatrix = TRUE)
#gset <- getGEO("GSE13159", GSEMatrix =TRUE, getGPL=FALSE)
my_gse <- if(is.list(gse)) gse[[3]] else gse
data_matrix <- exprs(my_gse)
data=t(data_matrix)

if (length(gse) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gse[[idx]]

meta_data <- pData(my_gse)

dim(meta_data)

colnames(meta_data)

head(meta_data[, ])
table(meta_data$`leukemia class:ch1`)
meta_data$Group <- ifelse(meta_data$`leukemia class:ch1` == "Non-leukemia and healthy bone marrow", 
                          "Control", 
                          "Disease")

meta_data$Group <- factor(meta_data$Group, levels = c("Control", "Disease"))

table(meta_data$Group)

"213261_at" %in% colnames(data)
x.comp <- data[, "213261_at"]
plot(x.comp~meta_data$Group)

##selected features
target_probes <- c("213261_at", "221979_at", "209586_s_at",  
                   "219966_x_at", "212560_at", "200736_s_at", "212549_at", 
                   "210183_x_at", "218999_at", "203288_at", "78383_at", 
                   "222307_at", "213506_at")

present_probes <- target_probes[target_probes %in% colnames(data)]
missing_probes <- target_probes[!target_probes %in% colnames(data)]

sub_data <- data[, present_probes, drop = FALSE]


library(RColorBrewer)
my_colors <- brewer.pal(n = 9, name = "RdYlBu")
col_palette <- colorRampPalette(my_colors)(200) 
par(mfrow=c(1,2), mar = 0.1+c(4,4,2,1))
cor1=cor(sub_data)
corrplot::corrplot(cor1, 
                   method = "color",     
                   col = col_palette,    
                   tl.col = "black",
                   mar = c(4, 0, 0, 0),
                   type = "lower", 
                   
)

fit=glm(y3~sub_data,family=binomial())
summary(fit)
p.fit=fit$fitted.values
library(pROC)
roc_obj=roc(y3,as.vector(p.fit),plot=TRUE,print.auc=TRUE,main='Ordinary LOgistic Regression',ci=T,
          col="lightpink3",lty=1,size=3,percent=T,direction="<")
auc=auc(y3,as.vector(p.fit))

library(glmnet)

cvfit1=cv.glmnet(sub_data[,],y3, family="binomial", alpha=0.5,nfolds=10 )
lambda_best1 <- cvfit1$lambda.min
b1=coef(cvfit1,s=lambda_best1)
preds1 <- predict(cvfit1, newx=as.matrix(sub_data[,]), s=lambda_best1, type="response")
roc_obj1=roc(y3,as.vector(preds1),plot=TRUE,print.auc=TRUE,main='Elastic Net Logistic Regression',ci=T
          ,col="lightpink3",lty=1,size=3,percent=T,direction="<")
auc1=auc(y3,as.vector(preds1))



