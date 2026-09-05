#### Data collection and processing ####
#GEO
library(tidyverse)
chooseBioCmirror()
BiocManager::install('GEOquery')
library(GEOquery)
gset = getGEO('GSE104948', destdir=".", AnnotGPL = F, getGPL = F)
#gset = getGEO('GSE40568', destdir=".", AnnotGPL = F, getGPL = F)
gset[[1]]
pdata <- pData(gset[[1]])
table(pdata$characteristics_ch1.1)
library(stringr)
row_tom1 <- filter(pdata,characteristics_ch1.1 == "diagnosis: ANCA Associated Vasculitis")
row_tom2 <- filter(pdata,characteristics_ch1.1== "")
data <- rbind(row_tom1, row_tom2)
exp <- exprs(gset[[1]])
comgene <- intersect(colnames(exp),rownames(data))
expr <- exp[,comgene]
boxplot(expr,outline=FALSE, notch=T,col=group_list, las=2)
dev.off()
###数据校正
library(limma)
exp=normalizeBetweenArrays(expr)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=3)
range(exp)
exp1 <- log2(exp+1)
range(exp1)
exp=exp1
dev.off()
qx <- as.numeric(quantile(exp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { exp[which(exp <= 0)] <- NaN
exprSet <- log2(exp)
print("log2 transform finished")}else{print("log2 transform not needed")}
index = gset[[1]]@annotation
#Import the GPL and add comments
gpl <- na.omit(gpl)
exp <- as.data.frame(exp)
rownames(gpl) <- gpl$ID
gpl1 <- gpl[,-1]
comname<-intersect(rownames(exp),rownames(gpl))
exp <- exp[comname,]
gpl1 <- gpl[comname,]
exp1 <- cbind(gpl1,exp)
exp1 <- exp1[!duplicated(exp1$Symbol),]
rownames(exp1) <- exp1$Symbol
exp1 <- exp1[,-(1:4)]
write.table(exp1, file = "exp1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#IgG4-RD: Repeat the steps above

#SVA
library(sva)
exp1 <- read.table("Aexp1.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp1=as.matrix(exp1)
exp2 <- read.table("Iexp1.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
exp2=as.matrix(exp2)

library(sva)
library(limma)
library(dplyr)

common_genes <- intersect(rownames(exp1), rownames(exp2))
exp_ANCA <- exp1[common_genes, ]
exp_IgG4 <- exp2[common_genes, ]

merged_expr <- cbind(exp_ANCA, exp_IgG4)

# ANCA
pheno_ANCA <- data.frame(
  sample = colnames(exp_ANCA),
  disease = rep(c("Control", "ANCA"), times = c(18, 22)),
  tissue = rep("Kidney", 40),
  dataset = rep("GSE104948", 40),
  stringsAsFactors = FALSE
)

# IgG4
pheno_IgG4 <- data.frame(
  sample = colnames(exp_IgG4),
  disease = rep(c("IgG4", "Control"), times = c(5, 3)),
  tissue = rep("Salivary", 8),
  dataset = rep("GSE40568", 8),
  stringsAsFactors = FALSE
)

# Merge
pheno <- rbind(pheno_ANCA, pheno_IgG4)
table(pheno$disease, pheno$tissue)

pheno$disease <- factor(pheno$disease, levels = c("Control", "ANCA", "IgG4"))
mod <- model.matrix(~ disease, data = pheno)
mod0 <- model.matrix(~ 1, data = pheno)
merged_expr <- as.matrix(merged_expr)
svobj <- sva(merged_expr, mod, mod0)
print(paste("Estimated number of proxy variables", svobj$n.sv))

corrected_expr <- removeBatchEffect(
  merged_expr,
  covariates = svobj$sv,
  design = mod
)

# ANCA
anca_cols <- which(pheno$dataset == "GSE104948")
# IgG4
igg4_cols <- which(pheno$dataset == "GSE40568")

# Split
expr_ANCA_corrected <- corrected_expr[, anca_cols]
expr_IgG4_corrected <- corrected_expr[, igg4_cols]

# pheno
pheno_ANCA <- pheno[anca_cols, ]
pheno_IgG4 <- pheno[igg4_cols, ]

write.table(expr_ANCA_corrected, file = "expr_ANCA_corrected.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(expr_IgG4_corrected, file = "expr_IgG4_corrected.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(pheno_ANCA, file = "pheno_ANCA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
write.table(pheno_IgG4, file = "pheno_IgG4.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#Visualization
library(ggplot2)
library(patchwork)
library(ggpubr)

pca_before <- prcomp(t(merged_expr), scale. = TRUE)
pca_after  <- prcomp(t(corrected_expr), scale. = TRUE)
pca_before_df <- as.data.frame(pca_before$x[, 1:2])
pca_after_df  <- as.data.frame(pca_after$x[, 1:2])

pca_before_df$disease <- pheno$disease
pca_before_df$tissue  <- pheno$tissue
pca_after_df$disease  <- pheno$disease
pca_after_df$tissue   <- pheno$tissue

var_before <- round(summary(pca_before)$importance[2, 1:2] * 100, 2)
var_after  <- round(summary(pca_after)$importance[2, 1:2] * 100, 2)

# PCA Plot Before Calibration
p1 <- ggplot(pca_before_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = disease, shape = tissue), size = 4, stroke = 1) +
  scale_color_manual(values = c("ANCA" = "red", "IgG4" = "orange", "Control" = "blue")) +
  scale_shape_manual(values = c("Kidney" = 16, "Salivary" = 17)) +
  labs(
    title = "PCA Plot Before Calibration",
    x = paste0("PC1 (", var_before[1], "%)"),
    y = paste0("PC2 (", var_before[2], "%)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# PCA Plot After Calibration
p2 <- ggplot(pca_after_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = disease, shape = tissue), size = 4, stroke = 1) +
  scale_color_manual(values = c("ANCA" = "red", "IgG4" = "orange", "Control" = "blue")) +
  scale_shape_manual(values = c("Kidney" = 16, "Salivary" = 17)) +
  labs(
    title = "PCA Plot After Calibration)",
    x = paste0("PC1 (", var_after[1], "%)"),
    y = paste0("PC2 (", var_after[2], "%)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )
print(p1)
print(p2)


####Identification of DEGs####
# ANCA
pheno_ANCA$disease <- factor(pheno_ANCA$disease,levels = c("Control", "ANCA"))
design_ANCA <- model.matrix(~ disease, data = pheno_ANCA)
fit_ANCA <- lmFit(expr_ANCA_corrected, design_ANCA)
fit_ANCA <- eBayes(fit_ANCA)
DEGs_ANCA <- topTable(fit_ANCA, coef  = "diseaseANCA", number = Inf)
write.table(DEGs_ANCA, file = "DEGs_ANCA.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# IgG4
design_IgG4 <- model.matrix(~ disease, data = pheno_IgG4)
fit_IgG4 <- lmFit(expr_IgG4_corrected, design_IgG4)
fit_IgG4 <- eBayes(fit_IgG4)
DEGs_IgG4 <- topTable(fit_IgG4, coef = "diseaseIgG4", number = Inf)
write.table(DEGs_IgG4, file = "DEGs_IgG4.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
# Identify Up- and Down-Regulated Genes(ANCA)
degANCA=DEGs_ANCA
logFC=1
P.Value = 0.05
k1 = (degANCA$P.Value< P.Value)&(degANCA$logFC < -logFC)
k2 = (degANCA$P.Value < P.Value)&(degANCA$logFC > logFC)
degANCA$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(degANCA$change)
write.table(degANCA, file = "degANCA_change.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# # Identify Up- and Down-Regulated Genes(IgG4)
degIgG4=DEGs_IgG4
logFC=1
P.Value = 0.05
k1 = (degIgG4$P.Value< P.Value)&(degIgG4$logFC < -logFC)
k2 = (degIgG4$P.Value < P.Value)&(degIgG4$logFC > logFC)
degIgG4$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(degIgG4$change)

#ANCA heatmaps
library(pheatmap)
choose_gene=head(rownames(DEGs_ANCA),50)
choose_matrix=expr_ANCA_corrected[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix)
p1 <-pheatmap(choose_matrix,
              color = colorRampPalette(c("navy", "white", "red"))(100),
              breaks = seq(-3, 3, length.out = 101))
print(p1)

#IgG4 heatmaps
library(pheatmap)
choose_gene=head(rownames(DEGs_IgG4),50)
choose_matrix=expr_IgG4_corrected[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
p2 <-pheatmap(choose_matrix,
              color = colorRampPalette(c("navy", "white", "red"))(100),
              breaks = seq(-3, 3, length.out = 101))
print(p2)


#volcano plots
library(ggplot2)
this_tile <- paste0('Volcano plot of AAV',round(logFC,3),
                    '\nThe number of up gene is ',nrow(degANCA[degANCA$change =='up',]) ,
                    '\nThe number of down gene is ',nrow(degANCA[degANCA$change =='down',])
)
g1 = ggplot(data=degANCA,
            aes(x=logFC, y=-log10(P.Value),
                color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_bw(base_size = 10) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile ) +
  theme(plot.title = element_text(size=10,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))
print(g1)


this_tile <- paste0('Volcano plot of IgG4-RD',round(logFC,3),
                    '\nThe number of up gene is ',nrow(degIgG4[degIgG4$change =='up',]) ,
                    '\nThe number of down gene is ',nrow(degIgG4[degIgG4$change =='down',])
)
g1 = ggplot(data=degIgG4,
            aes(x=logFC, y=-log10(P.Value),
                color=change)) +
  geom_point(alpha=0.4, size=1.75) +
  theme_bw(base_size = 10) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle( this_tile ) +
  theme(plot.title = element_text(size=10,hjust = 0.5))+
  scale_colour_manual(values = c('blue','black','red'))
print(g1)

####Functional enrichment analysis####
library(tidyverse)
library("BiocManager")
library(org.Hs.eg.db)
library(clusterProfiler)
DEGs_ANCA <- read.table("degANCA_change.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
DEGs_ANCA <- DEGs_ANCA %>% rownames_to_column("Gene")
genelist <- bitr(DEGs_ANCA$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEGs_ANCA <- inner_join(DEGs_ANCA,genelist,by=c("Gene"="SYMBOL"))
###GO
ego <- enrichGO(gene = DEGs_ANCA$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "ALL",
                pAdjustMethod = "BH",
                minGSSize = 10,
                pvalueCutoff =0.05,
                qvalueCutoff =0.05,
                readable = TRUE)
ego_res <- ego@result
write.table(ego_res, file = "ego_res.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# Visualization
barplot(ego, drop = TRUE, showCategory =10,split="ONTOLOGY") +
  facet_grid(ONTOLOGY~., scale='free')

dev.off()
#KEGG
kk <- enrichKEGG(gene         = DEGs_ANCA$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.1,
                 qvalueCutoff =0.1)
kk_res <- kk@result
write.table(kk_res, file = "kk_res.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
# Visualization
barplot(kk, showCategory = 10,color = "pvalue")
dev.off()


####Lasso and Boruta####
#Lasso
library("glmnet")

gene=PPIMCC15
hubgenes=c(gene$V1)
hubgenes_expression<-expr_ANCA_corrected[match(hubgenes,rownames (expr_ANCA_corrected)),]
hubgenes_expression = na.omit(hubgenes_expression)
x=as.matrix(hubgenes_expression[,c(1:ncol(hubgenes_expression))])
set.seed(1)
design_ANCA=as.data.frame(design_ANCA)

y=data.matrix(design_ANCA$diseaseANCA)
x=t(x)

#LASSO regression with 10‑fold cross‑validation
fit=glmnet(x,y,family = "binomial",maxit = 1000)
plot(fit,xvar="lambda",label = TRUE)
cvfit = cv.glmnet(x,y,family="binomial",maxit = 1000)
plot(cvfit)
dev.off()


coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
#Output Results
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
geneCoef

#boruta
library(tidyverse)
library(Boruta)
library(caret)
colnames(group) <- c('sample', 'group')
gene=PPIMCC15
colnames(gene) <- "symbol"
dat <- expr_ANCA_corrected[rownames(expr_ANCA_corrected) %in% gene$symbol, ] %>%  t() %>%  as.data.frame()
dat$sample <- rownames(dat)
dat <- merge(dat, group, var = "sample")
dat <- column_to_rownames(dat, var = "sample") %>% as.data.frame()
table(dat$group)
dat$group <- factor(dat$group, levels = c('ANCA', 'control'))
#Analysis
set.seed(123)
boruta.train <- Boruta(group~., data = dat, doTrace = 2, maxRuns = 500)
final.boruta <- TentativeRoughFix(boruta.train)
boruta_gene <- data.frame(symbol = getSelectedAttributes(final.boruta, withTentative = T))
write.csv(boruta_gene, file = 'boruta_gene.csv')
# Boruta
{
  plot(final.boruta, xlab = "", xaxt = "n")
  lz<-lapply(1:ncol(final.boruta$ImpHistory),function(i)
    final.boruta$ImpHistory[is.finite(final.boruta$ImpHistory[,i]),i])
  names(lz) <- colnames(final.boruta$ImpHistory)
  Labels <- sort(sapply(lz,median))
  axis(side = 1,las=2,labels = names(Labels),
       at = 1:ncol(final.boruta$ImpHistory),
       cex.axis = 0.8)
}



#cox
library(ggrepel)
library(ggplot2)
library(ggpubr)
colors=c('#4169E1','#FF0000','#A136A1')
x=t(x)
write.table(design, file = "design.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- group %>% rownames_to_column("sample")

#ITGAM
dataITGAM= data.frame(t(x["ITGAM",group$sample]))#
dataITGAM=t(dataITGAM)
dataITGAM=data.frame(dataITGAM)

group_name =unique(group$group)
my_comparison =as.list(as.data.frame(combn(group_name,2)))
dataITGAM$Group = factor(group$group)
# Visualization(ITGAM)
ggplot(dataITGAM,aes(x =Group,y=dataITGAM,color=Group))+
  scale_color_manual(values=colors[1:2])+
  geom_boxplot(width=0.5)+
  geom_point(position = position_jitter(0.2),size=1)+
  stat_summary(fun=median,geom='crossbar',size=1,width=0.3)+
  stat_compare_means(comparisons=my_comparison,size=13,#p.format displays only the p-value and not the test method; p.signif displays the significance level symbol.
                     method='wilcox.test',exact=FALSE,aes(label=..p.signif..))+
  coord_cartesian(ylim = c(5,10))+
  theme_bw()+xlab("")+ylab("ITGAMGene(log2CPM)")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(size=20,colour="black",face="bold"),
        axis.title=element_text(size=20,colour="black",face= "bold"),
        legend.text=element_text(size=20,face="bold"),
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x= element_text(angle=0,hjust = 0.5))
dev.off()
write.csv(dataITGAM,file = "dataITGAM.csv",row.names = T)

#CXCL1
dataCXCL1= data.frame(t(x["CXCL1",group$sample]))#
dataCXCL1=t(dataCXCL1)
dataCXCL1=data.frame(dataCXCL1)
group_name =unique(group$group)
my_comparison =as.list(as.data.frame(combn(group_name,2)))
dataCXCL1$Group = factor(group$group)
# Visualization(CXCL1)
ggplot(dataCXCL1,aes(x =Group,y=dataCXCL1,color=Group))+
  scale_color_manual(values=colors[1:2])+
  geom_boxplot(width=0.5)+
  geom_point(position = position_jitter(0.2),size=1)+
  stat_summary(fun=median,geom='crossbar',size=1,width=0.3)+
  stat_compare_means(comparisons=my_comparison,size=13,
                     method='wilcox.test',exact=FALSE,aes(label=..p.signif..))+
  coord_cartesian(ylim = c(5,9))+
  theme_bw()+xlab("")+ylab("CXCL1Gene(log2CPM)")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(size=20,colour="black",face="bold"),
        axis.title=element_text(size=20,colour="black",face= "bold"),
        legend.text=element_text(size=20,face="bold"),
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x= element_text(angle=0,hjust = 0.5))
dev.off()
write.csv(dataCXCL1,file = "dataCXCL1.csv",row.names = T)

####ROC-Bootstrap####
hubgenes=c("ITGAM","CXCL1")
hubgenes_expression<-expr_ANCA_corrected[match(hubgenes,rownames (expr_ANCA_corrected)),]
hubgenes_expression=as.matrix(hubgenes_expression)
library(pROC)
library(pROC)
library(caret)
set.seed(100)

# Data
n <- ncol(hubgenes_expression)
gene_names <- rownames(hubgenes_expression)
expr_t <- t(hubgenes_expression)

set.seed(100)
B <- 200
boot_auc_itgam <- c()
boot_auc_cxcl1 <- c()
y_num <- as.numeric(as.character(y))
for (b in 1:B) {
  boot_idx <- sample(1:n, size = n, replace = TRUE)
  oob_idx <- setdiff(1:n, unique(boot_idx))
  if (length(oob_idx) < 3 || length(unique(y_num[oob_idx])) < 2) next
  # ITGAM
  train_df <- data.frame(y = y_num[boot_idx], gene = expr_t[boot_idx, "ITGAM"])
  if (length(unique(train_df$y)) > 1) {
    fit <- glm(y ~ gene, data = train_df, family = binomial, maxit = 1000)
    pred <- predict(fit, newdata = data.frame(gene = expr_t[oob_idx, "ITGAM"]), type = "response")
    if (length(pred) == length(oob_idx)) {
      boot_auc_itgam <- c(boot_auc_itgam, roc(y_num[oob_idx], pred, quiet = TRUE)$auc)
    }
  }
  # CXCL1
  train_df <- data.frame(y = y_num[boot_idx], gene = expr_t[boot_idx, "CXCL1"])
  if (length(unique(train_df$y)) > 1) {
    fit <- glm(y ~ gene, data = train_df, family = binomial, maxit = 1000)
    pred <- predict(fit, newdata = data.frame(gene = expr_t[oob_idx, "CXCL1"]), type = "response")
    if (length(pred) == length(oob_idx)) {
      boot_auc_cxcl1 <- c(boot_auc_cxcl1, roc(y_num[oob_idx], pred, quiet = TRUE)$auc)
    }
  }
}

cat("\n========== Bootstrap Final Results (200 resamples)==========\n")
cat(sprintf("ITGAM successfully performed %d calculations, with an average AUC: %.3f (95%% CI: %.3f - %.3f)\n",
            length(boot_auc_itgam),
            mean(boot_auc_itgam),
            quantile(boot_auc_itgam, 0.025),
            quantile(boot_auc_itgam, 0.975)))
cat(sprintf("CXCL1 successfully performed %d calculations, with an average AUC: %.3f (95%% CI: %.3f - %.3f)\n",
            length(boot_auc_cxcl1),
            mean(boot_auc_cxcl1),
            quantile(boot_auc_cxcl1, 0.025),
            quantile(boot_auc_cxcl1, 0.975)))
#ROC
library(pROC)
set.seed(100)
# Bootstrap
B <- 200
n <- nrow(expr_t)

all_y_oob_itgam <- c()
all_pred_oob_itgam <- c()
all_y_oob_cxcl1 <- c()
all_pred_oob_cxcl1 <- c()

for (b in 1:B) {
  boot_idx <- sample(1:n, size = n, replace = TRUE)
  oob_idx <- setdiff(1:n, unique(boot_idx))
  if (length(oob_idx) < 3 || length(unique(y_num[oob_idx])) < 2) next
  #ITGAM
  train_df <- data.frame(y = y_num[boot_idx], gene = expr_t[boot_idx, "ITGAM"])
  fit <- glm(y ~ gene, data = train_df, family = binomial, maxit = 1000)
  pred <- predict(fit, newdata = data.frame(gene = expr_t[oob_idx, "ITGAM"]), type = "response")
  all_y_oob_itgam <- c(all_y_oob_itgam, y_num[oob_idx])
  all_pred_oob_itgam <- c(all_pred_oob_itgam, pred)
  # CXCL1
  train_df <- data.frame(y = y_num[boot_idx], gene = expr_t[boot_idx, "CXCL1"])
  fit <- glm(y ~ gene, data = train_df, family = binomial, maxit = 1000)
  pred <- predict(fit, newdata = data.frame(gene = expr_t[oob_idx, "CXCL1"]), type = "response")
  all_y_oob_cxcl1 <- c(all_y_oob_cxcl1, y_num[oob_idx])
  all_pred_oob_cxcl1 <- c(all_pred_oob_cxcl1, pred)
}

roc_boot_itgam <- roc(all_y_oob_itgam, all_pred_oob_itgam)
roc_boot_cxcl1 <- roc(all_y_oob_cxcl1, all_pred_oob_cxcl1)

# AUC
auc_itgam_boot <- round(auc(roc_boot_itgam), 3)
auc_cxcl1_boot <- round(auc(roc_boot_cxcl1), 3)


plot(roc_boot_itgam, col = "red", legacy.axes = TRUE,
     main = "Bootstrap-Corrected ROC Curves")
plot(roc_boot_cxcl1, add = TRUE, col = "blue")
legend("bottomright",
       legend = c("ITGAM (AUC:0.921 95%CI:0.909 0.921 0.933)", "CXCL1 (AUC:0.937 95%CI:0.927 0.937 0.947)"),
       col = c("red", "blue"), lwd = 2)
#AUC and 95% IC
round(auc(roc_boot_itgam), 3)
round(ci(roc_boot_itgam), 3)
round(auc(roc_boot_cxcl1), 3)
round(ci(roc_boot_cxcl1), 3)
legend("bottomright",
       legend = c("ITGAM (AUC:0.921 95%CI:0.909 0.921 0.933)", "CXCL1 (AUC:0.937 95%CI:0.927 0.937 0.947)"),
       col = c("red", "blue"), lwd = 2)

####GSEA####
library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(tidyverse)
library(pacman)
library(clusterProfiler)
library(enrichplot)
ITGAM_exp <- expr_ANCA_corrected["ITGAM", ]
case_samples <- pheno_ANCA$sample[pheno_ANCA$disease == "ANCA"]
expr_case <- expr_ANCA_corrected[, case_samples]
ITGAM_exp_case <- expr_case["ITGAM", ]
ITGAM_exp_case <- as.numeric(ITGAM_exp_case)
group_ITGAM <- ifelse(ITGAM_exp_case > median(ITGAM_exp_case), "High", "Low")
group_ITGAM <- factor(group_ITGAM, levels = c("Low", "High"))

#Identification of DEGs
design <- model.matrix(~ 0 + group_ITGAM)
colnames(design) <- levels(group_ITGAM)
fit <- lmFit(expr_case, design)
contrast_matrix <- makeContrasts(High - Low, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

res_ITGAM <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")
geneList <- res_ITGAM$t#t-test
names(geneList) <- rownames(res_ITGAM)
geneList <- geneList[!is.na(names(geneList))]
geneList <- sort(geneList, decreasing = TRUE)

library(org.Hs.eg.db)
gene_symbols <- names(geneList)
entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_symbols, column = "ENTREZID",
                     keytype = "SYMBOL", multiVals = "first")

valid_idx <- !is.na(entrez_ids)
geneList_entrez <- geneList[valid_idx]
names(geneList_entrez) <- entrez_ids[valid_idx]

geneList_entrez <- geneList_entrez[!duplicated(names(geneList_entrez))]


set.seed(123)
gsea_results <- gseKEGG(
  geneList = geneList_entrez,
  organism = "hsa",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.25,
  eps = 0,
  verbose = FALSE
)


significant_pathways <- gsea_results[gsea_results$p.adjust < 0.25, ]
if (nrow(significant_pathways) > 10) {
  significant_pathways <- significant_pathways[1:10, ]  # top10
}
# Visualization
ridge_p <- ridgeplot(gsea_results, showCategory = 10, fill = "p.adjust") +
  theme_minimal(base_size = 12) +
  labs(title = "GSEA of ITGAM High vs Low Expression in ANCA",
       x = "Running Enrichment Score", y = "Pathways") +
  theme(axis.text.y = element_text(size = 10),
        legend.position = "right")


ggsave("Figure_GSEA_Ridgeplot_ITGAM.tiff", plot = ridge_p,
       width = 8, height = 10, dpi = 300, compression = "lzw")

print(significant_pathways[, c("Description", "NES", "p.adjust")])
write.table(significant_pathways, file = "significant_pathwaysITGAM_IGG4.txt",sep = "\t",row.names = T,col.names = NA,quote = F)



####Immune cell infiltration####
#cibersort
library(e1071)
library(parallel)
library(preprocessCore)
library(tidyverse)
source("CIBERSORT.R")
sig_matrix <- "LM22.txt"
mixture_file = 'expr_IgG4_corrected.txt'

res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
save(res_cibersort,file = "res_cibersort.Rdata")
load("res_cibersort.Rdata")
res_cibersort <- res_cibersort[,1:22]
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]
write.table(ciber.res, file = "CIBERSORT-Results.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

# Visualization
library(RColorBrewer)
comname<-intersect(rownames(ciber.res),rownames(group))
ciber.res <- ciber.res[comname,]
group <- group[comname,]
ciber.res1 <- cbind(ciber.res,group)
ciber.res1 <- as.data.frame(ciber.res1)
group_vec <- ciber.res1$group

n_types <- ncol(ciber.res1)
mycol <- ggplot2::alpha(hcl.colors(ncol(ciber.res), palette = "Dynamic"), 1)
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7)
par(bty = "o",
    mgp = c(5, 2.5, 0),
    mar = c(5.1, 6.0, 4.1, 14),
    tcl = -0.25,
    las = 1,
    xpd = FALSE)
bp <- barplot(as.matrix(t(ciber.res1)),
              border = NA,
              names.arg = rep("", nrow(ciber.res1)),
              yaxt = "n",
              ylab = "Relative percentage",
              col = mycol,
              main = "Immune Cell Abundance (CIBERSORT)",
              cex.main = 1.5)
axis(side = 2,
     at = c(0, 0.2, 0.4, 0.6, 0.8, 1),
     labels = c("0%", "20%", "40%", "60%", "80%", "100%"))
split_indices <- split(seq_along(bp), ciber.res1$group)
group_midpoints <- sapply(split_indices, function(idx) mean(bp[idx]))
group_colors <- c("black", "black")
names(group_colors) <- names(group_midpoints)
mapply(function(name, mid, col) {
  mtext(name, side = 1, line = 1, at = mid, cex = 1.2, font = 2, col = col)
}, names(group_midpoints), group_midpoints, group_colors[names(group_midpoints)])

legend(par("usr")[2] - 0,
       par("usr")[4],
       legend = colnames(ciber.res1),
       xpd = TRUE,
       fill = mycol,
       cex = 0.8,
       border = NA,
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()
# Visualization
library(tidyverse)
a <- read.table("CIBERSORT-Results.txt", sep = "\t",row.names = 1,check.names = F,header = T)
colnames(group) <- c('sample', 'group')
a$group <- group$group[match(rownames(a), group$sample)]
b <- group
a <- a %>% rownames_to_column("sample")
library(ggsci)
library(tidyr)
library(ggpubr)

b <- gather(a,key=CIBERSORT,value = Proportion,-c(group,sample))

ggboxplot(b, x = "CIBERSORT", y = "Proportion",
          fill = "group", palette = "lancet") +
  scale_fill_manual(values = c("#DF1C26", "navy")) +
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test",
    label = "p.signif",
    label.size = 5,  # 显著性标记（*、ns）字号，默认约3.88，调大至5
    symnum.args = list(
      cutpoints = c(0, 0.001, 0.01, 0.05, 0.075, 1),
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  theme(
    text = element_text(size = 14),          # 基础字号从10增至14
    axis.title = element_text(size = 14),    # 轴标题（可选，默认继承text）
    axis.text = element_text(size = 13),     # 轴刻度数字（略小于标题）
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13), # 保留角度
    legend.title = element_text(size = 14),  # 图例标题
    legend.text = element_text(size = 13),   # 图例条目
    plot.title = element_text(size = 16, face = "bold") # 若有主标题则调大
  )
dev.off()


#ssgsea
library(ggsci)
library(tidyr)
library(ggpubr)
library(tidyverse)
library(data.table)
library(GSEABase)
library(GSVA)
library(tibble)
library(ggpubr)
library(ggcorrplot)
library(ggplot2)
cellreports <- read.table("geneSet.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
list<-as.data.frame(t(cellreports))
list<-as.list(list)#list

expr_IgG4_corrected=read.table("expr_IgG4_corrected.txt",header=T,sep="\t",row.names=1)

params <- gsvaParam(as.matrix(expr_IgG4_corrected),list ,
                    minSize = 1,
                    maxSize = Inf,
                    kcdf = "Gaussian",
                    tau = 1,
                    maxDiff = TRUE,
                    absRanking = FALSE)
ssgsea <- gsva(params,
               verbose = TRUE,
               BPPARAM = BiocParallel::SerialParam(progressbar = TRUE))
library(pheatmap)
library(limma)
library(dplyr)
library(vegan)
library(ggplot2)
library(stringr)
expr_IgG4_corrected <- expr_IgG4_corrected %>% t() %>% as.data.frame()#t转换
load("design-group_list.Rda")
write.table(design, file = "design.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
group=read.table("group.txt",header=T,sep="\t",row.names=1)
annotation_col = data.frame(group_list)
rownames(annotation_col)=rownames(design)
colnames(annotation_col)[1]<-"group"
ssgsea.1 <-decostand(ssgsea,"standardize",MARGIN =1)
#heatmap
p=pheatmap(
  ssgsea.1,
  border_color = NA,
  cluster_rows = T,cluster_cols = F,
  color = colorRampPalette(colors = c("blue","white","tomato"))(100),
  labels_row = NULL,
  annotation_col = annotation_col,
  clustering_method = "ward.D2",
  fontsize_col = 3,
  cutree_cols = 2,
  show_rownames = T,
  show_colnames = F,
)
y=t(ssgsea.1)
data <- cbind(y,group)
data <-data[,c(29,30,1:28)]
data=pivot_longer(data=data,
                  cols = 3:30,
                  names_to = "celltype",
                  values_to = "proportion")
ggboxplot(data = data,
          x = "celltype",
          y = "proportion",
          combine = FALSE,
          merge = FALSE,
          color = "black",
          fill = "group",
          palette = c("navy","#DF1C26"),
          title = NULL,
          xlab = "ssGSEA",
          ylab = "Expression",
          bxp.errorbar = FALSE,
          bxp.errorbar.width = 0.2,
          facet.by = NULL,
          panel.labs = NULL,
          short.panel.labs = TRUE,
          linetype = "solid",
          size = NULL,
          width = NULL,
          notch = FALSE,
          outlier.shape = 20,
          select = NULL,
          remove = NULL,
          order = NULL,
          error.plot = "pointrange",
          label = NULL,
          font.label = list(size = 12, color = "black"),
          label.select = NULL,
          repel = TRUE,
          label.rectangle = TRUE,
          ggtheme = theme_pubr()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))

#MCPcounter
library(curl)
library(devtools)
library(MCPcounter)
genes <- data.table::fread("genes.txt",data.table = F)
probesets <- data.table::fread("probesets.txt",data.table = F,header = F)
exp1=read.table("expr_ANCA_corrected.txt",header=T,sep="\t",row.names=1)

results<- MCPcounter.estimate(exp1,
                              featuresType="HUGO_symbols",
                              probesets=probesets,
                              genes=genes)
write.table(results,file="MCPCounter_scoreIgG4.txt",sep="\t",row.names =T,quote=F)

# heatmap
sample_names <- colnames(results)
group_vec <- group$group[match(sample_names, group$sample)]
group_colors <- ifelse(group_vec == "control", "grey","red")
heatmap(as.matrix(results),
        col = colorRampPalette(c("blue", "white", "red"))(100),
        ColSideColors = group_colors,
        scale = "row",
        cexRow = 0.6,
        cexCol = 0.8
)
library(reshape)
results = as.data.frame(results)
results$Cell = rownames(results)
group=read.table("group.txt",header=T,sep="\t",row.names=1)
colnames(group) <- c('sample', 'group')

cell_long <- merge(group, melt(results, id = "Cell"),
                   by.x = "sample", by.y = "variable")
library(ggplot2)
library(tibble)
library(ggpubr)
library(pheatmap)
ggplot(cell_long, aes(Cell, value, fill = group)) +
  geom_boxplot(outlier.shape = 21, color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 80, vjust = 0.5)) +
  stat_compare_means(aes(group = group, label = ..p.signif..),
                     method = "kruskal.test",
                     symnum.args = list(
                       cutpoints = c(0, 0.001, 0.01, 0.05, 0.075, Inf),
                       symbols = c("****", "***", "**", "*", "ns")
                     ))

mcp_score <- as.data.frame(t(results))
head(mcp_score)
overlap_genes <- c("ITGAM","CXCL1")
gene_expression <- t(exp1[overlap_genes, ])
combined_data <- merge(gene_expression, mcp_score, by = "row.names")
cell_cols <- colnames(mcp_score)
combined_data[cell_cols] <- lapply(combined_data[cell_cols], as.numeric)

#Correlation Calculation
cor_results <- data.frame()
for (gene in overlap_genes) {
  for (cell in cell_cols) {
    cor_test <- cor.test(
      combined_data[[gene]],
      combined_data[[cell]],
      method = "spearman"
    )
    cor_results <- rbind(cor_results, data.frame(
      Gene = gene,
      CellType = cell,
      Rho = cor_test$estimate,
      Pvalue = cor_test$p.value
    ))
  }
}

# ================== 新增 FDR 校正 ================== #
# 使用 Benjamini-Hochberg 方法控制错误发现率
cor_results$FDR <- p.adjust(cor_results$Pvalue, method = "BH")
# ================================================== #

# 筛选显著相关（基于 FDR < 0.05）
significant_cor <- cor_results[cor_results$FDR < 0.05, ]

# 后续绘图：建议使用 FDR 标记显著性
library(ggplot2)
library(dplyr)

# 修改显著性标记列（使用 FDR）
cor_results <- cor_results %>%
  mutate(
    neg_log_p = -log10(Pvalue),  # 气泡大小仍用原始 P 值
    Sig = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01  ~ "**",
      FDR < 0.05  ~ "*",
      TRUE        ~ ""
    )
  )

# 气泡图（使用 FDR 标记）
p_bubble <- ggplot(cor_results, aes(x = CellType, y = Gene)) +
  geom_point(aes(color = Rho, size = neg_log_p), shape = 18) +
  geom_text(aes(label = Sig), color = "black", size = 4, vjust = 0.8) +
  scale_color_gradient2(
    low = "navy",
    mid = "white",
    high = "#DF1C26",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  scale_size_continuous(range = c(3, 10), name = "-log10(Pvalue)") +
  labs(
    title = "Spearman Correlation: Genes vs Immune Cells",
    x = "Cell Types",
    y = "Genes",
    color = "Spearman Rho"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )
print(p_bubble)

# 热图（使用 FDR 标记）
p_heatmap <- ggplot(cor_results, aes(x = CellType, y = Gene, fill = Rho)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = paste0(round(Rho, 2), "\n", Sig)),
            color = "black", size = 3.5) +
  scale_fill_gradient2(
    low = "#3498DB",
    mid = "white",
    high = "#E74C3C",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  labs(title = "Spearman Correlation Heatmap", x = "Cell Types", y = "Genes") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
print(p_heatmap)



####GSE108109####
setwd("GSE108109")
library(tidyverse)
chooseBioCmirror()
BiocManager::install('GEOquery')
library(GEOquery)
gset = getGEO('GSE108109', destdir=".", AnnotGPL = F, getGPL = F)

gset[[1]]
pdata <- pData(gset[[1]])
table(pdata$source_name_ch1)
library(stringr)
row_tom1 <- filter(pdata,source_name_ch1 == "Kidney Biopsy from Human ANCA Associated Vasculitis")
row_tom2 <- filter(pdata,source_name_ch1== "Kidney Biopsy from Human Living donor")
data <- rbind(row_tom1, row_tom2)


group_list <- ifelse(str_detect(data$source_name_ch1, "Kidney Biopsy from Human ANCA Associated Vasculitis"), "ANCA",
                     "control")
group_list = factor(group_list,
                    levels = c("control","ANCA"))

exp <- exprs(gset[[1]])
comgene <- intersect(colnames(exp),rownames(data))
expr <- exp[,comgene]
boxplot(expr,outline=FALSE, notch=T,col=group_list, las=2)
dev.off()

library(limma)
exp=normalizeBetweenArrays(expr)
boxplot(exp,outline=FALSE, notch=T,col=group_list, las=3)
range(exp)
exp1 <- log2(exp+1)
range(exp1)
exp=exp1
dev.off()
qx <- as.numeric(quantile(exp, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { exp[which(exp <= 0)] <- NaN
exprSet <- log2(exp)
print("log2 transform finished")}else{print("log2 transform not needed")}
index = gset[[1]]@annotation
#Import the GPL and add comments
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
colnames(gpl)
x1=gpl$ENTREZ_GENE_ID
x1=as.character (x1)
x2=AnnotationDbi::select(org.Hs.eg.db , keys=x1,
                         columns = c("ENTREZID" , "SYMBOL"),keytype = "ENTREZID")
gpl1<- gpl[, c(1, 2)]
gpl1$gene=x2$SYMBOL
gpl <- gpl1

gpl <- na.omit(gpl)
exp <- as.data.frame(exp)
rownames(gpl) <- gpl$ID
gpl1 <- gpl[,-1]
comname<-intersect(rownames(exp),rownames(gpl1))
exp <- exp[comname,]
gpl1 <- gpl1[comname,]
exp1 <- cbind(gpl1,exp)
exp1 <- exp1[!duplicated(exp1$gene),]
rownames(exp1) <- exp1$gene
exp1 <- exp1[,-(1:2)]
write.table(exp1, file = "exp1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
#Corresponding Groups
identical(colnames(exp1), rownames(data))
data <- data[colnames(exp1), ]
group_list <- ifelse(str_detect(data$source_name_ch1, "Kidney Biopsy from Human ANCA Associated Vasculitis"), "ANCA",
                     "control")
group_list = factor(group_list,
                    levels = c("control","ANCA"))

#Identification of DEGs
library(limma)
design=model.matrix(~0+factor(group_list))
rownames(design)=rownames(data)
colnames(design)=levels(factor(group_list))

contrast.matrix=makeContrasts(ANCA-control,levels = design)
fit=lmFit(exp1,design)
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)
deg=topTable(fit,coef=1,number = Inf)
nrDEG = na.omit(deg)
head(nrDEG)
write.table(nrDEG, file = "deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
nrDEG <- read.table("deg_all.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
deg=nrDEG
#Identify Up- and Down-Regulated Genes
logFC=1
P.Value = 0.05
k1 = (deg$P.Value< P.Value)&(deg$logFC < -logFC)
k2 = (deg$P.Value < P.Value)&(deg$logFC > logFC)
deg$change = ifelse(k1,"down",ifelse(k2,"up","stable"))
table(deg$change)
write.table(deg, file = "deg_change.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

hubgenes=c("ITGAM","CXCL1")
hubgenes_expression<-exp1[match(hubgenes,rownames (exp1)),]
hubgenes_expression = na.omit(hubgenes_expression)
x=as.matrix(hubgenes_expression[,c(1:ncol(hubgenes_expression))])
#cox
library(ggrepel)
library(ggplot2)
library(ggpubr)
colors=c('#4169E1','#FF0000','#A136A1')
write.table(design, file = "design.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
group <- read.table("group.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
group <- group %>% rownames_to_column("sample")

#ITGAM
dataITGAM= data.frame(t(x["ITGAM",group$sample]))#
dataITGAM=t(dataITGAM)
dataITGAM=data.frame(dataITGAM)

group_name =unique(group$group)
my_comparison =as.list(as.data.frame(combn(group_name,2)))
dataITGAM$Group = factor(group$group)
# Visualization(ITGAM)
ggplot(dataITGAM,aes(x =Group,y=dataITGAM,color=Group))+
  scale_color_manual(values=colors[1:2])+
  geom_boxplot(width=0.5)+
  geom_point(position = position_jitter(0.2),size=1)+
  stat_summary(fun=median,geom='crossbar',size=1,width=0.3)+
  stat_compare_means(comparisons=my_comparison,size=13,#p.format displays only the p-value and not the test method; p.signif displays the significance level symbol.
                     method='wilcox.test',exact=FALSE,aes(label=..p.signif..))+
  coord_cartesian(ylim = c(0,5))+
  theme_bw()+xlab("")+ylab("ITGAMGene(log2CPM)")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(size=20,colour="black",face="bold"),
        axis.title=element_text(size=20,colour="black",face= "bold"),
        legend.text=element_text(size=20,face="bold"),
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x= element_text(angle=0,hjust = 0.5))
dev.off()
write.csv(dataITGAM,file = "dataITGAM.csv",row.names = T)

#CXCL1
dataCXCL1= data.frame(t(x["CXCL1",group$sample]))#
dataCXCL1=t(dataCXCL1)
dataCXCL1=data.frame(dataCXCL1)
group_name =unique(group$group)
my_comparison =as.list(as.data.frame(combn(group_name,2)))
dataCXCL1$Group = factor(group$group)
# Visualization(CXCL1)
ggplot(dataCXCL1,aes(x =Group,y=dataCXCL1,color=Group))+
  scale_color_manual(values=colors[1:2])+
  geom_boxplot(width=0.5)+
  geom_point(position = position_jitter(0.2),size=1)+
  stat_summary(fun=median,geom='crossbar',size=1,width=0.3)+
  stat_compare_means(comparisons=my_comparison,size=13,
                     method='wilcox.test',exact=FALSE,aes(label=..p.signif..))+
  coord_cartesian(ylim = c(1,9))+
  theme_bw()+xlab("")+ylab("CXCL1Gene(log2CPM)")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(size=20,colour="black",face="bold"),
        axis.title=element_text(size=20,colour="black",face= "bold"),
        legend.text=element_text(size=20,face="bold"),
        legend.title=element_blank(),
        legend.key=element_blank(),
        axis.text.x= element_text(angle=0,hjust = 0.5))
dev.off()
write.csv(dataCXCL1,file = "dataCXCL1.csv",row.names = T)

#ROC
hubgenes=c("ITGAM","CXCL1")
hubgenes_expression<-exp1[match(hubgenes,rownames (exp1)),]
hubgenes_expression=as.matrix(hubgenes_expression)
library(pROC)
roc1<- roc(group_list, hubgenes_expression[1,])
roc2<- roc(group_list, hubgenes_expression[2,])
plot(roc1, col = "red", legacy.axes = TRUE,
     main = "External validation ROC curves")
plot(roc2, add=TRUE, col="blue")
round(auc(roc1),3)##AUC
round(ci(roc1),3)##95%CI
round(auc(roc2),3)##AUC
round(ci(roc2),3)##95%CI
legend("bottomright",
       legend = c("ITGAM (AUC:0.844 95%CI:0.673 0.844 1.000)", "CXCL1 (AUC:0.622 95%CI:0.341 0.622 0.903)"),
       col = c("red", "blue"), lwd = 2)
dev.off()
