library(Hmisc)
library(StepReg)
library(pROC)
#######READ DATA########
##You can download data directly from GEO 
##https://ftp.ncbi.nlm.nih.gov/geo/series/GSE9nnn/GSE9476/matrix/ 
## Or from Supporting file S1
####1-From GEO###########
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

library(GEOquery)
getOption("download.file.method.GEOquery")
options("download.file.method.GEOquery" = "libcurl") 
gse<- getGEO(filename = "C:///Downloads/GSE9476_series_matrix.txt.gz", getGPL = F, GSEMatrix = TRUE)
my_gse <- if(is.list(gse)) gse[[3]] else gse
data_matrix <- exprs(my_gse)
data=t(data_matrix)
##2- From supporting file S1 
data <- read.delim("C:///Downloads/S1 File.txt", header=FALSE)
data <- read.delim(file.choose(), header=FALSE)

##1- data prparation
x=as.matrix(t(data[,-1]))
sum(is.na(x))
colnames(x)=as.vector(data[,1])
y=factor(c(rep(0,18),rep(1,26), rep(0,20))) 
data.lukemia=data.frame(y,x)

## subset selcetion in data
y3=y[19:64]
x3=x[c(19:64),]

##2 find the complete separation
min0=max0=0
for ( i in 1:ncol(x3) ) {
  if (min(x3[y3==0,i]) > max(x3[y3==1,i])) {min0[i]=i}
  if (max(x3[y3==0,i])<min(x3[y3==1,i])) {max0[i]=i}
}

table(min0);table(max0)
sum(min0!=0,na.rm=T);sum(max0!=0,na.rm=T)
id.comp=sum(min0,na.rm=T)
x.comp=as.vector(x3[,id.comp])
names(x.comp)=colnames(x3)[id.comp]
y3=factor(y3, labels = c("Normal","AML"))
table(y3)





df <- data.frame(x =x.comp, y = y3)
library(ggplot2)

p1<-ggplot(df, aes(x = x, color = y)) +
  geom_vline(xintercept = 6.36, 
             linetype = "dashed",  
             color = "black",      
             linewidth = 0.8) +    
  stat_ecdf(linewidth = 1) + 
  theme_minimal() +
  labs(title = "Empirical Distribution Function", x = "213261_at",
       y = "F(x)",
       color = "Groups") +
  scale_color_manual(values = c("Normal" = "#1f77b4", "AML" = "#ff7f0e")) 
p2 <- ggplot(df, aes(x = y, y = x, fill = y)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 6.36, 
             linetype = "dashed",  
             color = "black",      
             linewidth = 0.8) + 
  theme_minimal() +
  labs(title = "Boxplot",
     x="", y = "213261_at") +
  scale_fill_manual(values = c("Normal" = "#1f77b4", "AML" = "#ff7f0e")) +
  guides(fill = "none") 
p3=ggplot(df, aes(x = x, fill = y)) +
  geom_density(alpha = 0.4) + 
  theme_minimal() +
  labs(title = "Density Plot",
       x = "213261_at",
       y = "f(x)",
       fill = "Groups") +
  scale_fill_manual(values = c("Normal" = "#1f77b4", "AML" = "#ff7f0e")) 

library(patchwork)
combined_p <- p2+p1/p3
print(combined_p)
ggsave("E://papers/MLS/PLOS one/Figs/Fig2.png", plot = combined_p, width = 8, height = 6, dpi = 600)
##dropping complete separation from features

x3 <- x3[, -id.comp]

##3-Mann-Whitney Test
p=0
for(i in 1:ncol(x3)){
  p[i]=wilcox.test(x3[,i]~y3,exact=FALSE)$p.value
}

##4 Feature Selection
m<-cbind(1:ncol(x3),p)
msort=m[order(m[,2],decreasing=FALSE),]	
p.fdr=0
for(i in 1:nrow(m)){
  p.fdr[i]=msort[i,2]*nrow(m)/i#.05#.015
}
p.fdr[1:30]
k.FDR=cumsum(p.fdr<.0005);k.FDR
k.fdr=0
for (i in 1:nrow(m)){
  k.fdr[i]= (ifelse (k.FDR[i+1]==k.FDR[i], i ,NA) )
}
k.FDR=min(k.fdr,na.rm=T);k.FDR
xk.FDR=x3[,msort[1:k.FDR,1]]
dim(xk.FDR)
cors1=cor(xk.FDR)
l=0
for(j in 1:(k.FDR-1)){
  l[j]=(1-abs(cors1[j, j+1]))/(1+abs(cors1[j, j+1]))
}
l=c(1,l)
p.fdr3=0
for(i in 1:k.FDR){
  p.fdr3[i]=msort[i,2]*nrow(m)/sum(l[1:i]) #.05#.015
}
k.FDR3=cumsum(p.fdr3<.0005)
k.fdr3=0
for (i in 1:nrow(m)){
  k.fdr3[i]= (ifelse (k.FDR3[i+1]==k.FDR3[i], i ,NA) )
}
k1=min(k.fdr3,na.rm=T);k1

xk1=x3[,msort[1:k1,1]]
df=data.frame(x=xk1,y=y3)
par(mfrow = c(4,3), bty="o", mgp = c(2,.8,0), mar = 0.1+c(2,4,2,1))
p=list()
for(i in 1:k1){
 
    temp_df <- data.frame(
      y = df$y,          
      val = xk1[, i]     
      )
    
       current_y_label <- colnames(xk1)[i]
    
       p[[i]] <- ggplot(temp_df, aes(x = y, y = val, fill = y)) +
         geom_boxplot(alpha = 0.7) +
         theme_minimal() +
         theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8), 
    panel.grid.major = element_line(colour = "gray90"),
    panel.grid.minor = element_blank()
  ) +
         labs(x = "", 
              y = current_y_label) + 
         scale_fill_manual(values = c("Normal" = "#1f77b4", "AML" = "#ff7f0e")) +
         guides(fill = "none")
}
library(patchwork)

combined_features <- wrap_plots(p, ncol = 4, nrow = 3)
print(combined_features)
ggsave("E://papers/MLS/PLOS one/Figs/Fig3.png", plot = combined_features, width = 8, height = 6, dpi = 600)

##5- Hierarchical Clustering of Features

xk1=x3[,msort[1:k1,1]]
corxk1=cor(cbind("213261_at"=x.comp,xk1))
library(RColorBrewer)
my_colors <- brewer.pal(n = 9, name = "RdYlBu")
col_palette <- colorRampPalette(my_colors)(200) 
par(mfrow=c(1,2), mar = 0.1+c(4,4,2,1))
corrplot::corrplot(corxk1, 
                   method = "color",     
                   col = col_palette,    
                   tl.col = "black",
                   mar = c(4, 0, 0, 0),
                   type = "lower", 
                   
)
mtext("Pearson Correlation", side = 1, line = 2, cex = 1, font = 2)


library(dendextend)
gg <- varclus(xk1, similarity="pearson", type="data.matrix", trans="square")
n_clus=cutree(gg$hclust,h=0.75^2)
gv <- as.dendrogram(gg$hclust)
gv <- color_branches(gv, k = 3) 
gv <- set(gv, par = list(lwd = 5))
par(mar = 0.1+c(6,4,2,1),bty="7",col.axis="black",col.lab="black")
plot(gv, horiz = FALSE, axes = FALSE,  ylim = c(0, .8),lwd = 5,font.lab=2)
abline(h=0.75^2,col="red3",lty=3)
rh <- seq(0, 1, by=0.1)  # re-label x-axis re:similarity not distance
axis(2, at=1 - rh, labels=format(rh))
box()
mtext("1 - R Square", side = 2, line = 2, cex = 1, font = 2)

clus=as.numeric(n_clus)
length(n_clus)
table(n_clus)
xkmem=rbind(t(clus),xk1)
xkmem1=xkmem[-1,xkmem[1,]==1]
xkmem2=xkmem[-1,xkmem[1,]==2]
xkmem3=xkmem[-1,xkmem[1,]==3]

colnames=c(colnames(xkmem1),colnames(xkmem2),colnames(xkmem3))
#write.csv(colnames,"H:/papers/MLS/luk/Data& Codes/colnames.csv")

##6- Logitic Regression by Elastic Net Penalties
library(glmnet)


preds0=1-(x.comp/max(x.comp))
auc0=auc(y3,as.vector(x.comp))

cvfit1=cv.glmnet(xkmem1,y3, family="binomial", alpha=0.5,nfolds=10 )
lambda_best1 <- cvfit1$lambda.min
preds1 <- predict(cvfit1, newx=as.matrix(xkmem1), s=lambda_best1, type="response")
#pROC::roc(y3,as.vector(preds1),plot=TRUE,print.auc=TRUE,main='CLUSTER=1',col="lightpink3",lty=1,size=3,percent=T,direction="<")
auc1=auc(y3,as.vector(preds1))
b1=coef(cvfit1,s=lambda_best1)


cvfit2=cv.glmnet(xkmem2,y3, family="binomial", alpha=0.5,nfolds=10 )
lambda_best2 <- cvfit2$lambda.min
preds2 <- predict(cvfit2, newx=as.matrix(xkmem2), s=lambda_best2, type="response")
#pROC::roc(y3,as.vector(preds2),plot=TRUE,print.auc=TRUE,main='CLUSTER=2',col="lightpink3",lty=1,size=3,percent=T,direction="<")
auc2=auc(y3,as.vector(preds2))
b2=coef(cvfit2,s=lambda_best2)


cvfit3=cv.glmnet(xkmem3,y3, family="binomial", alpha=0.5,nfolds=10 )
lambda_best3 <- cvfit3$lambda.min
preds3 <- predict(cvfit3, newx=as.matrix(xkmem3), s=lambda_best3, type="response")
auc3=auc(y3,as.vector(preds3))
b3=coef(cvfit3,s=lambda_best3)


table1=rbind(b1,b2,b3)
t1=table1[,1]
t2=row.names(table1)
#write.csv(cbind(t2,t1),"H:/papers/MLS/luk/Data& Codes/table1.csv")


##7- Ensemble of Models

enspred=rowMeans(cbind(preds0,preds1,preds2,preds3))
auc.ens=auc(y3,as.vector(enspred))

##8- Bootstrap
t.sys0=Sys.time()
xkmem.sim=rbind(t(clus),xk1)
y.sim=y[c(19:64)]
t1=table(y.sim)[2];t0=table(y.sim)[1]
nsim=1000
auc1=auc2=auc3=auc.ens=auc.comp=0
b1=matrix(0,nrow=ncol(xkmem1)+1,ncol=nsim)
b2=matrix(0,nrow=ncol(xkmem2)+1,ncol=nsim)
b3=matrix(0,nrow=ncol(xkmem3)+1,ncol=nsim)
preds0=0
for(bs in 1:nsim){
  id.sample1=sample(1:t1,t1,replace =T)
  id.sample0=sample((t1+1):(t1+t0),t0,replace =T)
  id.sample=c(id.sample1,id.sample0)
  y3=y.sim[id.sample]
  xkmem=xkmem.sim[c(1,id.sample),]
  dim(xkmem)
  xkmem1=xkmem[-1,xkmem[1,]==1]
  xkmem2=xkmem[-1,xkmem[1,]==2]
  xkmem3=xkmem[-1,xkmem[1,]==3]
  
  ## logitic regression by elasticnet approach
  
  cvfit1=cv.glmnet(xkmem1,y3, family="binomial", alpha=0.5,nfolds=10 )
  lambda_best1 <- cvfit1$lambda.min
  preds1 <- predict(cvfit1, newx=as.matrix(xkmem1), s=lambda_best1, type="response")
  #pROC::roc(y3,as.vector(preds1),plot=TRUE,print.auc=TRUE,main='CLUSTER=1',col="blue3")
  auc1[bs]=auc(y3,as.vector(preds1))
  b1.1=coef(cvfit1,s=lambda_best1)
  b1[,bs]=b1.1[,1]
  
  
  cvfit2=cv.glmnet(xkmem2,y3, family="binomial", alpha=0.5,nfolds=10 )
  lambda_best2 <- cvfit2$lambda.min
  preds2 <- predict(cvfit2, newx=as.matrix(xkmem2), s=lambda_best2, type="response")
  #pROC::roc(y3,as.vector(preds2),plot=TRUE,print.auc=TRUE,main='CLUSTER=2',col="blue3")
  auc2[bs]=auc(y3,as.vector(preds2))
  b2.1=coef(cvfit2,s=lambda_best2)
  b2[,bs]=b2.1[,1]
  
  
  cvfit3=cv.glmnet(xkmem3,y3, family="binomial", alpha=0.5,nfolds=10 )
  lambda_best3 <- cvfit3$lambda.min
  preds3 <- predict(cvfit3, newx=as.matrix(xkmem3), s=lambda_best3, type="response")
  #pROC::roc(y3,as.vector(preds3),plot=TRUE,print.auc=TRUE,main='CLUSTER=3',col="blue3")
  auc3[bs]=auc(y3,as.vector(preds3))
  b3.1=coef(cvfit3,s=lambda_best3)
  b3[,bs]=b3.1[,1]
  
  preds0=1-(x.comp/max(x.comp))
  auc.comp[bs]=auc(y3,as.vector(preds0))
  
  
  enspred=rowMeans(cbind(preds0,preds1,preds2,preds3))
  #pROC::roc(y3,as.vector(enspred),plot=TRUE,print.auc=TRUE,main='Ensemble',col="blue3")
  auc.ens[bs]=auc(y3,as.vector(enspred))
}

mean(auc1); mean(auc2); mean(auc3); mean(auc.comp);mean(auc.ens)

betas=t(rbind(b1[-1,],b2[-1,],b3[-1,]))
colnames(betas)=colnames
betai=se.betai=lcibeta=ucibeta=oddsi=lciodds=uciodds=0
for(i in 1:ncol(betas)){
  betai[i]=mean(betas[,i])
  se.betai[i]=sd(betas[,i])/sqrt(nrow(x3))
  lcibeta[i]=betai[i]+qnorm(0.00025)*se.betai[i]
  ucibeta[i]=betai[i]+qnorm(0.99975)*se.betai[i]
  oddsi[i]=exp(betai[i])
  lciodds[i]=exp(lcibeta[i])
  uciodds[i]=exp(ucibeta[i])
}
table3=round(cbind(betai,se.betai,lcibeta,ucibeta,oddsi,lciodds,uciodds),3)
me_genes <- data.frame(
  gene = colnames,
  mean = betai,
  lower = lcibeta,
  upper = ucibeta,
  Clusters = c("cluster 1","cluster 1","cluster 1","cluster 1","cluster 1","cluster 1","cluster 2",
               "cluster 2","cluster 2","cluster 2","cluster 3","cluster 3")
)
library(dplyr)
me_genes <- me_genes %>%
    arrange(mean) %>%
    mutate(gene = factor(gene, levels = gene))

IM.me <- ggplot(me_genes, aes(x = mean, y = gene)) +
  geom_errorbarh(aes(xmin = lower, xmax = upper, color = Clusters), height = 0.3, size = 0.8) +
  geom_point(aes(color = Clusters), size = 2.5) +
  
  scale_color_manual(values = c("cluster 1" = "forestgreen", "cluster 2" = "firebrick", "cluster 3" = "blue")) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.7) +
  theme_bw() + 
  labs(
    title = "Proposed Algorithm",
    x = "Cluster-Wise Elastic Net Logistic Regression Coefficient with CI",
    y = "Gene ID"
  ) +
    theme(
    panel.background = element_rect(fill = "White", color = NA),
    plot.background = element_rect(fill = "White", color = NA),
    panel.grid.major = element_line(color ="gray90", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = c(0.8, 0.2), 
    legend.background = element_rect(fill = "White", color = "black"), 
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
    
  )

print(IM.me)
ggsave("E://papers/MLS/PLOS one/Figs/Fig5.png", plot = IM.me, width = 8, height = 6, dpi = 600)

#########################################################################################
p_df=data.frame(preds0,preds1,preds2,preds3,enspred)
colnames(p_df)=c("Separation", "Cluster1", "Cluster2", "Cluster3", "Ensemble")
p_roc_ci <- list()
for(i in 1:5) {
  
  
  temp_df <- data.frame(
    outcome = df$y,
    predictor = p_df[, i]
  )
  
    roc_obj <- roc(temp_df$outcome, temp_df$predictor, quiet = TRUE)
    auc_val <- round(auc(roc_obj), 3)
    current_var_name <- colnames(p_df)[i] 

ci_auc <- ci.auc(roc_obj,conf.level=0.9995)
ci_text <- paste0("AUC = ", auc_val, ", CI [", round(ci_auc[1], 3), "-", round(ci_auc[3], 3), "]")

roc_data <- data.frame(
  Sensitivity = roc_obj$sensitivities,
  Specificity = roc_obj$specificities,
  FPR = 1 - roc_obj$specificities
)

p_roc_ci[[i]] <- ggplot(roc_data, aes(x = FPR, y = Sensitivity)) +
  geom_path(color = "#1f77b4", linewidth = 1) + 
  geom_abline(linetype = "dashed", color = "red", alpha = 0.5) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    subtitle =paste("ROC:", c(current_var_name, ci_text)), 
    x = ci_text,
    y = ""
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
}


##############################
aucs=c(auc.comp,auc1,auc2,auc3,auc.ens)
AUC_names= rep(c("C.Separation", "Cluster1", "Cluster2", "Cluster3", "Ensemble"), each=nsim)
df.auc=data.frame(aucs,group=AUC_names)
p6 <- ggplot(df.auc, aes(x = group, y = aucs, fill = group)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8), 
    panel.grid.major = element_line(colour = "gray90"),
    panel.grid.minor = element_blank()
  ) +
  labs(subtitle = "Boxplot of Bootstrap AUCs ",
       x="", y = "AUC") +
  scale_fill_manual(values = c("Complete Separation" = "#1f77b4","Cluster1"="#1f77b4",
                               
                               "Cluster2"="#1f77b4","Cluster3"="#1f77b4", "Ensemble" = "#ff7f0e")) +
  guides(fill = "none") 


combined_roc_ci <- wrap_plots(p_roc_ci, ncol = 2) + p6+
  plot_annotation(title = "ROC Curves ",
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5)))

print(combined_roc_ci)
ggsave("E://papers/MLS/PLOS one/Figs/Fig6.png", plot = combined_roc_ci, width = 8, height = 6, dpi = 600)

write.csv(table3,"E:/papers/MLS/luk/Data& Codes/table3.csv")
##################################################################################
##################################################################################

