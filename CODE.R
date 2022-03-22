###############   00.数据处理####
setwd("E:/JZH/00.data/GSE30153/")
tmp=read.table("GSE30153_series_matrix.id.txt",
               header=T,sep="\t",check.names = F)
tmp_uniq <- aggregate(.~id,tmp,max)###取平均值
write.table(tmp_uniq,file="GSE30153_series_matrix.id.uniq.txt",
            sep = '\t',row.names = F,quote = F)


setwd("E:/JZH/00.data/GSE65391/")
tmp=read.table("GSE65391_series_matrix.id.txt",
               header=T,sep="\t",check.names = F)
tmp_uniq <- aggregate(.~id,tmp,max)###取平均值
write.table(tmp_uniq,file="GSE23205_series_matrix.id.uniq.txt",
            sep = '\t',row.names = F,quote = F)

tmp1=read.table("GSE30153_series_matrix.id.uniq.group.txt",
                header=TRUE,check.names = FALSE, row.names = 1)
tmp2 <- data.frame(t(tmp1),check.names = F)
tmp3 <- tmp2[order(tmp2$group,decreasing = F),]
tmp4 <- data.frame(t(tmp3),check.names = F)
write.table(tmp4,file="GSE30153_series_matrix.id.uniq.group.sort.tmp.txt",sep = '\t',
            row.names = T,quote = F)



setwd("E:/JZH/00.data/GSE72509/")
tmp=read.table("GSE72509_SLE_RPKMs.id.txt",
               header=T,sep="\t",check.names = F)
tmp_uniq <- aggregate(.~id,tmp,max)###取平均值
write.table(tmp_uniq,file="GSE72509_SLE_RPKMs.id.uniq.txt",
            sep = '\t',row.names = F,quote = F)

tmp1=read.table("GSE72509_SLE_RPKMs.id.uniq.group.txt",
                header=TRUE,check.names = FALSE, row.names = 1)
tmp2 <- data.frame(t(tmp1),check.names = F)
tmp3 <- tmp2[order(tmp2$group,decreasing = F),]
tmp4 <- data.frame(t(tmp3),check.names = F)
write.table(tmp4,file="GSE72509_SLE_RPKMs.id.uniq.group.sort.tmp.txt",sep = '\t',
            row.names = T,quote = F)



#########   01.差异分析 ####
######### 01.1 GSE30153####
setwd("E:/JZH/01.DEG/GSE30153/")
library('gplots')
library('limma')
library(dplyr)
mRNA=read.table("../../00.data/GSE30153/GSE30153_series_matrix.id.uniq.txt",
                header=TRUE,row.names=1,check.names = FALSE)
par(mfrow=c(1,2))
#boxplot(data.frame(mRNA),col="blue")    ####画箱式图，比较数据分布情况，数据分布好，则不用进行log2转换
#mRNA <- log2(mRNA+1)
mRNA[1:5,1:5]
mRNA.group <- c(rep("normal_A",9),rep('tomor_N',17)) %>% factor(.,levels = c("normal_A","tomor_N"),ordered = F)
#group <- group[,1] #定义比较组，按照癌症和正常样品数#
mRNA.group <- model.matrix(~factor(mRNA.group))#把group设置成一个model matrix#
mRNA.fit <- lmFit(mRNA,mRNA.group)
mRNA.fit <- eBayes(mRNA.fit)
tempOutput = topTable(mRNA.fit,coef=2,n=Inf,adjust="BH")
nrDEG = na.omit(tempOutput)
mRNA.diff <- nrDEG
write.csv(mRNA.diff, "limmaOut.csv")
foldChange=1
padj=0.05
mRNA.diffSig = mRNA.diff[(mRNA.diff$P.Value < padj & (mRNA.diff$logFC>foldChange | mRNA.diff$logFC<(-foldChange))),]
write.csv(mRNA.diffSig, "diffSig.csv")
mRNA.diffUp = mRNA.diff[(mRNA.diff$P.Value < padj & (mRNA.diff$logFC>foldChange)),]
write.csv(mRNA.diffUp, "diffUp.csv")
mRNA.diffDown = mRNA.diff[(mRNA.diff$P.Value < padj & (mRNA.diff$logFC<(-foldChange))),]
write.csv(mRNA.diffDown, "diffDown.csv")
logFC <-mRNA.diff$logFC
deg.padj <- mRNA.diff$P.Value
data <- data.frame(logFC=logFC,padj=deg.padj)
data$mRNA.diffsig[(data$padj > 0.05|data$padj=="NA")|(data$logFC < foldChange)& data$logFC -foldChange] <- "no"
data$mRNA.diffsig[data$padj <= 0.05 & data$logFC >= foldChange] <- "up"
data$mRNA.diffsig[data$padj <= 0.05 & data$logFC <= -foldChange] <- "down"
x_lim <- max(logFC,-logFC)
library(ggplot2)
library(RColorBrewer)
pdf(file = "volcano.pdf",width=8,height=8)
theme_set(theme_bw())
p <- ggplot(data,aes(logFC,-1*log10(padj),color=mRNA.diffsig))+geom_point()+
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)
p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))
p <- p  +guides(colour = FALSE)
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
p
dev.off()
print(p)

##########   差异热图
library(RColorBrewer)
library(pheatmap)
gene.a <- read.table("heatmap.txt",header=T,sep="\t",check.names = FALSE,row.names = 1)
annotation_col=read.table("anno.txt",header=TRUE,
                          row.names=1,check.names = FALSE,sep = '\t')
gene.a <- log2(gene.a+1)
pheatmap(gene.a,clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         fontsize = 10,clustering_method = "ward.D",
         annotation_col =annotation_col,
         show_rownames = T, #不显示行名
         show_colnames = T,
         treeheight_row = 50,treeheight_col = 50,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),
         scale= "row",border_color = NA, cluster_cols = FALSE)


#########  01.2 GSE65391####
setwd("E:/JZH/01.DEG/GSE65391/")
library('gplots')
library('limma')
library(dplyr)
mRNA=read.table("../../00.data/GSE65391/GSE65391_series_matrix.id.uniq.group.sort.txt",
                header=TRUE,row.names=1,check.names = FALSE)
par(mfrow=c(1,2))
#boxplot(data.frame(mRNA),col="blue")    ####画箱式图，比较数据分布情况，数据分布好，则不用进行log2转换
#mRNA <- log2(mRNA+1)
mRNA[1:5,1:5]
mRNA.group <- c(rep("normal_A",48),rep('tomor_N',924)) %>% factor(.,levels = c("normal_A","tomor_N"),ordered = F)
#group <- group[,1] #定义比较组，按照癌症和正常样品数#
mRNA.group <- model.matrix(~factor(mRNA.group))#把group设置成一个model matrix#
mRNA.fit <- lmFit(mRNA,mRNA.group)
mRNA.fit <- eBayes(mRNA.fit)
tempOutput = topTable(mRNA.fit,coef=2,n=Inf,adjust="BH")
nrDEG = na.omit(tempOutput)
mRNA.diff <- nrDEG
write.csv(mRNA.diff, "limmaOut.csv")
foldChange=1
padj=0.05
mRNA.diffSig = mRNA.diff[(mRNA.diff$P.Value < padj & (mRNA.diff$logFC>foldChange | mRNA.diff$logFC<(-foldChange))),]
write.csv(mRNA.diffSig, "diffSig.csv")
mRNA.diffUp = mRNA.diff[(mRNA.diff$P.Value < padj & (mRNA.diff$logFC>foldChange)),]
write.csv(mRNA.diffUp, "diffUp.csv")
mRNA.diffDown = mRNA.diff[(mRNA.diff$P.Value < padj & (mRNA.diff$logFC<(-foldChange))),]
write.csv(mRNA.diffDown, "diffDown.csv")
logFC <-mRNA.diff$logFC
deg.padj <- mRNA.diff$P.Value
data <- data.frame(logFC=logFC,padj=deg.padj)
data$mRNA.diffsig[(data$padj > 0.05|data$padj=="NA")|(data$logFC < foldChange)& data$logFC -foldChange] <- "no"
data$mRNA.diffsig[data$padj <= 0.05 & data$logFC >= foldChange] <- "up"
data$mRNA.diffsig[data$padj <= 0.05 & data$logFC <= -foldChange] <- "down"
x_lim <- max(logFC,-logFC)
library(ggplot2)
library(RColorBrewer)
pdf(file = "volcano.pdf",width=8,height=8)
theme_set(theme_bw())
p <- ggplot(data,aes(logFC,-1*log10(padj),color=mRNA.diffsig))+geom_point()+
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)
p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))
p <- p  +guides(colour = FALSE)
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
p
dev.off()
print(p)

##########   差异热图
library(RColorBrewer)
library(pheatmap)
gene.a <- read.table("heatmap.txt",header=T,sep="\t",check.names = FALSE,row.names = 1)
annotation_col=read.table("anno.txt",header=TRUE,
                          row.names=1,check.names = FALSE,sep = '\t')

gene.a1 <- apply(gene.a,2,function(x){as.numeric(x)})
row.names(gene.a1) <- row.names(gene.a)
gene.a1 = as.matrix(gene.a1)
gene.a1[which(gene.a1 > quantile(gene.a1,0.8))] <-  quantile(gene.a1,0.8)
gene.a1[which(gene.a1 < quantile(gene.a1,0.2))] <-  quantile(gene.a1,0.2)
pheatmap(gene.a1,clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         fontsize = 10,clustering_method = "ward.D",
         annotation_col =annotation_col,
         show_rownames = F, #不显示行名
         show_colnames = F,
         treeheight_row = 50,treeheight_col = 50,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),
         scale= "row",border_color = NA, cluster_cols = FALSE)






#########  01.3 GSE72509####
setwd("E:/JZH/01.DEG/GSE72509/")
library('gplots')
library('limma')
library(dplyr)
mRNA=read.table("../../00.data/GSE72509/GSE72509_SLE_RPKMs.id.uniq.group.sort.txt",
                header=TRUE,row.names=1,check.names = FALSE)
par(mfrow=c(1,2))
#boxplot(data.frame(mRNA),col="blue")    ####画箱式图，比较数据分布情况，数据分布好，则不用进行log2转换
mRNA <- log2(mRNA+1)
mRNA[1:5,1:5]
mRNA.group <- c(rep("normal_A",18),rep('tomor_N',99)) %>% factor(.,levels = c("normal_A","tomor_N"),ordered = F)
#group <- group[,1] #定义比较组，按照癌症和正常样品数#
mRNA.group <- model.matrix(~factor(mRNA.group))#把group设置成一个model matrix#
mRNA.fit <- lmFit(mRNA,mRNA.group)
mRNA.fit <- eBayes(mRNA.fit)
tempOutput = topTable(mRNA.fit,coef=2,n=Inf,adjust="BH")
nrDEG = na.omit(tempOutput)
mRNA.diff <- nrDEG
write.csv(mRNA.diff, "limmaOut.csv")
foldChange=1
padj=0.05
mRNA.diffSig = mRNA.diff[(mRNA.diff$P.Value < padj & (mRNA.diff$logFC>foldChange | mRNA.diff$logFC<(-foldChange))),]
write.csv(mRNA.diffSig, "diffSig.csv")
mRNA.diffUp = mRNA.diff[(mRNA.diff$P.Value < padj & (mRNA.diff$logFC>foldChange)),]
write.csv(mRNA.diffUp, "diffUp.csv")
mRNA.diffDown = mRNA.diff[(mRNA.diff$P.Value < padj & (mRNA.diff$logFC<(-foldChange))),]
write.csv(mRNA.diffDown, "diffDown.csv")
logFC <-mRNA.diff$logFC
deg.padj <- mRNA.diff$P.Value
data <- data.frame(logFC=logFC,padj=deg.padj)
data$mRNA.diffsig[(data$padj > 0.05|data$padj=="NA")|(data$logFC < foldChange)& data$logFC -foldChange] <- "no"
data$mRNA.diffsig[data$padj <= 0.05 & data$logFC >= foldChange] <- "up"
data$mRNA.diffsig[data$padj <= 0.05 & data$logFC <= -foldChange] <- "down"
x_lim <- max(logFC,-logFC)
library(ggplot2)
library(RColorBrewer)
pdf(file = "volcano.pdf",width=8,height=8)
theme_set(theme_bw())
p <- ggplot(data,aes(logFC,-1*log10(padj),color=mRNA.diffsig))+geom_point()+
  labs(x="log2(FoldChange)",y="-log10(Pvalue)")
p <- p + scale_color_manual(values =c("#0072B5","grey","#BC3C28"))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)
p <- p +theme(panel.grid =element_blank())+
  theme(axis.line = element_line(size=0))
p <- p  +guides(colour = FALSE)
p <- p +theme(axis.text=element_text(size=20),axis.title=element_text(size=20))
p
dev.off()
print(p)

library(RColorBrewer)
library(pheatmap)
gene.a <- read.table("heatmap.txt",header=T,sep="\t",check.names = FALSE,row.names = 1)
annotation_col=read.table("anno.txt",header=TRUE,
                          row.names=1,check.names = FALSE,sep = '\t')

gene.a <- log2(gene.a+1)
pheatmap(gene.a,clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         fontsize = 10,clustering_method = "ward.D",
         annotation_col =annotation_col,
         show_rownames = F, #不显示行名
         show_colnames = F,
         treeheight_row = 50,treeheight_col = 50,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),
         scale= "row",border_color = NA, cluster_cols = FALSE)


############################     02.RRA 2个数据集 ####
#install.packages("RobustRankAggreg")
library("RobustRankAggreg")
setwd("E:/JZH/02.RRA/")
meta <- read.delim("up.0.5.input.txt",na = "")
meta <- as.list(meta)
for (i in 1: length(meta)){meta[[i]] <- meta[[i]][!is.na(meta[[i]])]
meta[[i]] <- as.character(meta[[i]])}
count <- c(16181, 29738, 26366)
rankmat <- rankMatrix(meta, N = count)
ranks<- aggregateRanks(rmat=rankmat)
ranks$adjustedPval <-apply(cbind(ranks$Score * max(count), 1), 1, min)
up.results <- ranks[ranks$adjustedPval < 0.05,]
write.table(up.results, "up.results.0.5.txt", sep = "\t")

meta <- read.delim("down.0.5.input.txt",na = "")
meta <- as.list(meta)
for (i in 1: length(meta)){meta[[i]] <- meta[[i]][!is.na(meta[[i]])]
meta[[i]] <- as.character(meta[[i]])}
count <- c(16181, 29738, 26366)
rankmat <- rankMatrix(meta, N = count)
ranks<- aggregateRanks(rmat=rankmat)
ranks$adjustedPval <-apply(cbind(ranks$Score * max(count), 1), 1, min)
down.results <- ranks[ranks$adjustedPval < 0.05,]
write.table(down.results, "down.results.0.5.txt", sep = "\t")

##########   02.差异热图 ####
library(RColorBrewer)
library(pheatmap)
gene.a <- read.table("logfc热图.txt",header=T,sep="\t",check.names = FALSE,row.names = 1)
pheatmap(gene.a,
         #clustering_distance_rows = "correlation",
         #clustering_distance_cols = "correlation",
         fontsize = 10,
         #clustering_method = "ward.D",
         #annotation_col =annotation_col,
         show_rownames = T, #不显示行名
         show_colnames = T,
         treeheight_row = 50,treeheight_col = 50,
         color = colorRampPalette(brewer.pal(9,"Reds"))(100),
         #scale= "row",
         border_color = NA, cluster_cols = FALSE)




##############    03.富集分析 ####
setwd("E:/JZH/03.富集/")
library(org.Hs.eg.db)
library(stringr)
library(BiocGenerics)
library(clusterProfiler)
require(DOSE)
diffSig <- read.table("input.txt",header=F,sep="\t",row.names=1)
DEG.gene_symbol = as.character(rownames(diffSig))
DEG.entrez_id = mapIds(x = org.Hs.eg.db,keys = DEG.gene_symbol,
                       keytype = "SYMBOL",column = "ENTREZID")
gene = bitr(DEG.gene_symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
DEG.entrez_id = na.omit(DEG.entrez_id)
#write.table(DEG.entrez_id,"DEG.entrez_id.txt",sep="\t")

ego <- enrichGO(
  gene  = gene$ENTREZID,
  keyType = "ENTREZID",
  OrgDb   = org.Hs.eg.db,
  ont     = "all",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE)

pdf(file = "GO.bar.pdf",width=8,height=4)
barplot(ego, drop = TRUE, showCategory = 10,split="ONTOLOGY") +
  scale_y_discrete(labels=function(x) str_wrap(x, width=80))+
  facet_grid(ONTOLOGY~., scale='free')
dev.off()
write.table(ego,file="GO.txt",sep = '\t',row.names = F,quote = F)

ekegg <- enrichKEGG(
  gene = gene$ENTREZID,
  keyType = "kegg",
  organism  = 'hsa',
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH"
)

pdf(file = "KEGG.bar.pdf",width=8,height=2)
#dotplot(ekegg, orderBy='pvalue',decreasing=FALSE, showCategory = 10)
barplot(ekegg, orderBy='pvalue',decreasing=FALSE, showCategory = 10)
#dotplot(ekegg, showCategory = 10)
dev.off()
kegg.gene.list <- unlist(lapply(ekegg$geneID, function(v)
  paste(gene$SYMBOL[match(strsplit(v, '/', fixed=TRUE)[[1]],
                          gene$ENTREZID)], collapse = '/')))
write.table(kegg.gene.list,file="kegg.gene.list.txt",sep = '\t',row.names = F,quote = F)
write.table(ekegg,file="KEGG.txt",sep = '\t',row.names = F,quote = F)




#########################    04.机器学习####
#########   4.1 LR ####
set.seed(438)
setwd("E:/JZH/05.机器学习/")
library("pROC")
library(ggplot2)
library(ggthemes)
library(verification)
data_FPKM_TCGA <- read.table("input.txt",sep="\t",header=T,row.names=1)
geneExp_temp <- read.table("input.group.txt",sep="\t",header=T,row.names=1)
test <- read.table("测试集input.txt",sep="\t",header=T,row.names=1)
testgroup <- read.table("测试集.group.txt",sep="\t",header=T,row.names=1)


factor_name_list<-NULL
factor_name_ROC<-NULL
Factor_comb_matrix<-NULL
factor_comb_matrix.t <- NULL

geneExp<-data_FPKM_TCGA
for (f_num in 1:dim(geneExp)[2]) {
  Factor_comb_matrix<-t(combn(colnames(geneExp),f_num))  ##提取因子组合矩阵
  Factor_comb_matrix.t<-t(combn(colnames(test),f_num))
  for (i in 1:dim(Factor_comb_matrix)[1] ){
    factor_expre_matrix<-geneExp_temp[,c(as.vector(Factor_comb_matrix[i,]),"group")] ##把具体因子组合基因表达矩阵提取出来
    factor_expre_matrix.t<-testgroup[,c(as.vector(Factor_comb_matrix.t[i,]),"group")]
    LR_model <-glm(as.factor(group)~.,family=binomial(link = "logit"),data = factor_expre_matrix)
    glm.probs <- predict(LR_model,factor_expre_matrix.t[-dim(factor_expre_matrix.t)],type="response")
    factor_expre_matrix.t$model_score<-glm.probs
    roc.list <-roc(group~model_score,plot=F,ci=T,print.auc=TRUE,levels=c("0","1"),data = factor_expre_matrix.t)
    factor_name_list<-append(factor_name_list,paste(c(as.vector(Factor_comb_matrix.t[i,])),collapse=";"))
    factor_name_ROC<-append(factor_name_ROC,round(roc.list$auc,3))
  }
}

Factor_comb_ROC_matrix<-data.frame(factor_name_list,factor_name_ROC)
write.csv(Factor_comb_ROC_matrix,file="LR.output.csv",row.names = F,quote = F)



##############   4.2 ANN####
set.seed(438)
setwd("E:/JZH/05.机器学习/")
library(neuralnet)
data_FPKM_TCGA <- read.table("input.txt",sep="\t",header=T,row.names=1)
geneExp_temp <- read.table("input.group.txt",sep="\t",header=T,row.names=1)
test <- read.table("测试集input.txt",sep="\t",header=T,row.names=1)
testgroup <- read.table("测试集.group.txt",sep="\t",header=T,row.names=1)

factor_name_list<-NULL
factor_name_ROC<-NULL
Factor_comb_matrix<-NULL
factor_comb_matrix.t <- NULL
error_id <- NULL
set.seed(438)
geneExp<-data_FPKM_TCGA
for (f_num in 7:dim(geneExp)[2]) {
  Factor_comb_matrix<-t(combn(colnames(geneExp),f_num))  ##提取因子组合矩阵
  Factor_comb_matrix.t<-t(combn(colnames(test),f_num))
  for (i in 121:dim(Factor_comb_matrix)[1] ){
  #for (i in 15:25){
    factor_expre_matrix<-geneExp_temp[,c(as.vector(Factor_comb_matrix[i,]),"group")] ##把具体因子组合基因表达矩阵提取出来
    factor_expre_matrix.t<-testgroup[,c(as.vector(Factor_comb_matrix.t[i,]),"group")]
    x.inv <- try(ANN_model <-neuralnet(formula = group~.,
                          data = factor_expre_matrix,
                          hidden=1,
                          threshold=0.01,stepmax = 1e+05))
    if('Error' %in% class(x.inv)){
      error_id <- append(error_id,c(f_num,i))
      next
    }else{
      ANN_model <- x.inv
    }
    ANN.probs <- predict(ANN_model,factor_expre_matrix.t[-dim(factor_expre_matrix.t)],type="class")
    factor_expre_matrix.t$model_score<-as.numeric(ANN.probs)
    roc.list <-roc(group~model_score,plot=F,ci=T,print.auc=TRUE,levels=c("0","1"),data = factor_expre_matrix.t)
    factor_name_list<-append(factor_name_list,paste(c(as.vector(Factor_comb_matrix.t[i,])),collapse=";"))
    factor_name_ROC<-append(factor_name_ROC,round(roc.list$auc,3))
    print(data.frame(factor_name_list,factor_name_ROC))
  }
}

Factor_comb_ROC_matrix<-data.frame(factor_name_list,factor_name_ROC)
write.csv(Factor_comb_ROC_matrix,file="ANN.output.csv",row.names = F,quote = F)



#########   4.3 SVM ####
set.seed(438)
setwd("E:/JZH/05.机器学习/")
library("pROC")
library(e1071)
data_FPKM_TCGA <- read.table("input.txt",sep="\t",header=T,row.names=1)
geneExp_temp <- read.table("input.group.txt",sep="\t",header=T,row.names=1)
test <- read.table("测试集input.txt",sep="\t",header=T,row.names=1)
testgroup <- read.table("测试集.group.txt",sep="\t",header=T,row.names=1)

factor_name_list<-NULL
factor_name_ROC<-NULL
Factor_comb_matrix<-NULL
factor_comb_matrix.t <- NULL

geneExp<-data_FPKM_TCGA
for (f_num in 1:dim(geneExp)[2]) {
  Factor_comb_matrix<-t(combn(colnames(geneExp),f_num))  ##提取因子组合矩阵
  Factor_comb_matrix.t<-t(combn(colnames(test),f_num))
  for (i in 1:dim(Factor_comb_matrix)[1] ){
    factor_expre_matrix<-geneExp_temp[,c(as.vector(Factor_comb_matrix[i,]),"group")] ##把具体因子组合基因表达矩阵提取出来
    factor_expre_matrix.t<-testgroup[,c(as.vector(Factor_comb_matrix.t[i,]),"group")]
    svm_model <-svm(group ~ .,
                    data=factor_expre_matrix,
                    method = "svmRadial",
                    trControl = factor_expre_matrix,
                    tuneLength = 8,
                    metric = "ROC")
    glm.probs <- predict(svm_model,newdata = factor_expre_matrix.t[-dim(factor_expre_matrix.t)],
                         type="response")
    factor_expre_matrix.t$model_score<-glm.probs
    roc.list <-roc(group~model_score,plot=F,ci=T,print.auc=TRUE,levels=c("0","1"),data = factor_expre_matrix.t)
    factor_name_list<-append(factor_name_list,paste(c(as.vector(Factor_comb_matrix.t[i,])),collapse=";"))
    factor_name_ROC<-append(factor_name_ROC,round(roc.list$auc,3))
  }
}

Factor_comb_ROC_matrix<-data.frame(factor_name_list,factor_name_ROC)
write.csv(Factor_comb_ROC_matrix,file="SVM.output.csv",row.names = F,quote = F)



#########   4.4 Adaboost ####
set.seed(438)
setwd("E:/JZH/05.机器学习/")
library("pROC")
library(adabag)
data_FPKM_TCGA <- read.table("input.txt",sep="\t",header=T,row.names=1)
geneExp_temp <- read.table("input.group.txt",sep="\t",header=T,row.names=1)
test <- read.table("测试集input.txt",sep="\t",header=T,row.names=1)
testgroup <- read.table("测试集.group.txt",sep="\t",header=T,row.names=1)

factor_name_list<-NULL
factor_name_ROC<-NULL
Factor_comb_matrix<-NULL
factor_comb_matrix.t <- NULL

geneExp<-data_FPKM_TCGA
for (f_num in 1:dim(geneExp)[2]) {
  Factor_comb_matrix<-t(combn(colnames(geneExp),f_num))  ##提取因子组合矩阵
  Factor_comb_matrix.t<-t(combn(colnames(test),f_num))
  for (i in 1:dim(Factor_comb_matrix)[1] ){
    factor_expre_matrix<-geneExp_temp[,c(as.vector(Factor_comb_matrix[i,]),"group")] ##把具体因子组合基因表达矩阵提取出来
    factor_expre_matrix.t<-testgroup[,c(as.vector(Factor_comb_matrix.t[i,]),"group")]
    factor_expre_matrix$group <- as.factor(factor_expre_matrix$group)
    adaboost <-boosting(group ~ .,
                        data=factor_expre_matrix,
                        boos=TRUE, mfinal=10)
    glm.probs <- predict(adaboost,newdata = factor_expre_matrix.t[-dim(factor_expre_matrix.t)])$class
    factor_expre_matrix.t$model_score<-as.numeric(glm.probs)
    roc.list <-roc(group~model_score,plot=F,ci=T,print.auc=TRUE,levels=c("0","1"),data = factor_expre_matrix.t)
    factor_name_list<-append(factor_name_list,paste(c(as.vector(Factor_comb_matrix.t[i,])),collapse=";"))
    factor_name_ROC<-append(factor_name_ROC,round(roc.list$auc,3))
    print(c(factor_name_list,factor_name_ROC))
  }
}

Factor_comb_ROC_matrix<-data.frame(factor_name_list,factor_name_ROC)
write.csv(Factor_comb_ROC_matrix,file="Adaboost.output.csv",row.names = F,quote = F)




#########   4.5 RF ####
set.seed(438)
setwd("E:/JZH/05.机器学习/")
library("pROC")
library(randomForest)
data_FPKM_TCGA <- read.table("input.txt",sep="\t",header=T,row.names=1)
geneExp_temp <- read.table("input.group.txt",sep="\t",header=T,row.names=1)
test <- read.table("测试集input.txt",sep="\t",header=T,row.names=1)
testgroup <- read.table("测试集.group.txt",sep="\t",header=T,row.names=1)

factor_name_list<-NULL
factor_name_ROC<-NULL
Factor_comb_matrix<-NULL
Factor_comb_matrix.t <- NULL

geneExp<-data_FPKM_TCGA
for (f_num in 1:dim(geneExp)[2]) {
  Factor_comb_matrix<-t(combn(colnames(geneExp),f_num))  ##提取因子组合矩阵
  Factor_comb_matrix.t<-t(combn(colnames(test),f_num))
  for (i in 1:dim(Factor_comb_matrix)[1] ){
    factor_expre_matrix<-geneExp_temp[,c(as.vector(Factor_comb_matrix[i,]),"group")] ##把具体因子组合基因表达矩阵提取出来
    factor_expre_matrix.t<-testgroup[,c(as.vector(Factor_comb_matrix.t[i,]),"group")]
    RF_model <- randomForest(group ~ .,
                             data = factor_expre_matrix,
                             ntree = 500, # 树的数目，例文为1000
                             nPerm = 50, # 扰动次数，一般为50
                             mtry = floor(sqrt(ncol(factor_expre_matrix)-1)),
                             proximity = T,
                             importance = T)
    glm.probs <- predict(RF_model, newdata = factor_expre_matrix.t[-dim(factor_expre_matrix.t)],group="response")
    factor_expre_matrix.t$model_score<-as.numeric(glm.probs)
    roc.list <-roc(group~model_score,plot=F,ci=T,print.auc=TRUE,levels=c("0","1"),data = factor_expre_matrix.t)
    factor_name_list<-append(factor_name_list,paste(c(as.vector(Factor_comb_matrix.t[i,])),collapse=";"))
    factor_name_ROC<-append(factor_name_ROC,round(roc.list$auc,3))
  }
}

Factor_comb_ROC_matrix<-data.frame(factor_name_list,factor_name_ROC)
write.csv(Factor_comb_ROC_matrix,file="RF.output.csv",row.names = F,quote = F)


################    05.相关性 ####
setwd("E:/JZH/06.相关性/")
library(corrplot)
library(psych)
rt <- read.table("input.txt",sep="\t",header=T,row.names=1)
Cor <- corr.test(rt,method = "pearson",adjust = "fdr")
M <- Cor$r
p.mat <- Cor$p
head(M)
col2 <- colorRampPalette(c("#053061","#2166AC","#4393C3", "#92C5DE","#D1E5F0", "#FFFFFF",
                           "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F"))

write.csv(p.mat,file="p值.csv",row.names = F,quote = F)
write.csv(M,file="相关性值.csv",row.names = F,quote = F)

#相关性热图
pdf('Fig6.corrplot.pdf',width = 9,height = 9)
corrplot(M, type="full", order="hclust",
         p.mat = p.mat, sig.level = 0.05,
         insig = "blank",col=col2(200))
dev.off()



###############   06.验证表达量 ####
###############   6.1.外部验证 GSE39088 ####
setwd("E:/JZH/07.验证/GSE39088/")
tmp=read.table("10gene.id.txt",
               header=T,sep="\t",check.names = F)
tmp_uniq <- aggregate(.~id,tmp,max)###取平均值
write.table(tmp_uniq,file="10gene.id.uniq.txt",
            sep = '\t',row.names = F,quote = F)

library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
C <- read.table('box.txt',sep = '\t',header = T,check.names = F)
#C[,2:ncol(C)] <- apply(C[,2:ncol(C)],2,function(x){log2(x+1)})
colnames(C)[1] <- 'Type'
C1 <- gather(C,gene,expr,2:ncol(C))
ggboxplot(C1, x= 'gene', y='expr',
          ylab = "gene", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),label = 'p.signif',label.x = 1.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))####  横坐标倾斜45度


###############    6.2.整合数据集 box####
setwd("E:/JZH/07.验证/")
C <- read.table('训练集boxinput.txt',sep = '\t',header = T,check.names = F)
#C[,2:ncol(C)] <- apply(C[,2:ncol(C)],2,function(x){log2(x+1)})
colnames(C)[1] <- 'Type'
C1 <- gather(C,gene,expr,2:ncol(C))
ggboxplot(C1, x= 'gene', y='expr',
          ylab = "gene", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),label = 'p.signif',label.x = 1.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))####  横坐标倾斜45度



#################### 07.roc####
##################   7.1 训练集 ####
setwd("E:/JZH/08.ROC/")
data=read.table("整合数据集input.txt",header=T,sep="\t",row.names=1)
library("pROC")
mycol <- brewer.pal(10,'Set3')
#单基因ROC，画在一张图上
#data <- data.frame(t(data),check.names = F)
rocdata <- data.frame(Sample = rownames(data),
                      exp = data$EIF2AK2,
                      Type = data$lable)

library(RColorBrewer)
pdf('整合数据集ROC.pdf',width = 8,height = 8)
x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=F, #绘制平滑曲线
              main="",
              #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
              col=mycol[1],#线的颜色
              lwd=2, #线的粗细
              legacy.axes=T)
j=1
auc.test <- paste0(colnames(data)[1],' AUC : ',format(as.numeric(x$auc),digits=3))
for (i in colnames(data[3:ncol(data)-1])){
  j=j+1
  rocdata <- data.frame(Sample = rownames(data),
                        exp = data[,i],
                        Type = data$lable)
  x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[j],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T,add=T)
  auc.test <- c(auc.test,
                paste0(i,' AUC : ',format(as.numeric(x$auc),digits=3)))
}
legend(0.6,0.5, auc.test,lwd=2,bty="n",col=mycol,cex = 1.3)
dev.off()


##################   7.2 模型验证集 ####
setwd("E:/JZH/08.ROC/")
data=read.table("模型验证集input.txt",header=T,sep="\t",row.names=1)
library("pROC")
mycol <- brewer.pal(10,'Set3')
#单基因ROC，画在一张图上
#data <- data.frame(t(data),check.names = F)
rocdata <- data.frame(Sample = rownames(data),
                      exp = data$EIF2AK2,
                      Type = data$lable)

library(RColorBrewer)
pdf('模型验证集ROC.pdf',width = 8,height = 8)
x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=F, #绘制平滑曲线
              main="",
              #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
              col=mycol[1],#线的颜色
              lwd=2, #线的粗细
              legacy.axes=T)
j=1
auc.test <- paste0(colnames(data)[1],' AUC : ',format(as.numeric(x$auc),digits=3))
for (i in colnames(data[3:ncol(data)-1])){
  j=j+1
  rocdata <- data.frame(Sample = rownames(data),
                        exp = data[,i],
                        Type = data$lable)
  x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[j],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T,add=T)
  auc.test <- c(auc.test,
                paste0(i,' AUC : ',format(as.numeric(x$auc),digits=3)))
}
legend(0.6,0.5, auc.test,lwd=2,bty="n",col=mycol,cex = 1.3)
dev.off()


##################   7.3 验证集 ####
data=read.table("验证集input.txt",check.names =F,header=T,sep="\t",row.names=1)
library("pROC")
#单基因ROC，画在一张图上
#data <- data.frame(t(data),check.names = F)
rocdata <- data.frame(Sample = rownames(data),
                      exp = data$EIF2AK2,
                      Type = data$lable)

library(RColorBrewer)
pdf('验证集ROC.pdf',width = 8,height = 8)
x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=F, #绘制平滑曲线
              main="",
              #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
              col=mycol[1],#线的颜色
              lwd=2, #线的粗细
              legacy.axes=T)
j=1
auc.test <- paste0(colnames(data)[1],' AUC : ',format(as.numeric(x$auc),digits=3))
for (i in colnames(data[3:ncol(data)-1])){
  j=j+1
  rocdata <- data.frame(Sample = rownames(data),
                        exp = data[,i],
                        Type = data$lable)
  x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, #绘制平滑曲线
                main="",
                #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
                col=mycol[j],#线的颜色
                lwd=2, #线的粗细
                legacy.axes=T,add=T)
  auc.test <- c(auc.test,
                paste0(i,' AUC : ',format(as.numeric(x$auc),digits=3)))
}
legend(0.6,0.5, auc.test,lwd=2,bty="n",col=mycol,cex = 1.3)
dev.off()




###################    09.网络####
###################    10.预测药物####


#######################   ！！！！！！！！！！！11.xcell  免疫浸润！！！！！！！！！！！！！！！！！####
library(RColorBrewer)
library(pheatmap)
setwd("E:/JZH/11.免疫浸润xcell/")
gene.a <- read.table("heatmap.txt",header=T,sep="\t",check.names = FALSE,row.names = 1)
annotation_col=read.table("../00.data/GSE65391/anno.txt",header=TRUE,
                          row.names=1,check.names = FALSE,sep = '\t')

pheatmap(gene.a,clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         fontsize = 10,clustering_method = "ward.D",
         annotation_col = annotation_col,
         show_rownames = T, #不显示行名
        show_colnames = F,
         treeheight_row = 50,treeheight_col = 50,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),
         border_color = NA, cluster_cols = FALSE)


library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
C <- read.table('boxdata.txt',sep = '\t',header = T,check.names = F)
#C[,2:ncol(C)] <- apply(C[,2:ncol(C)],2,function(x){log2(x+1)})
colnames(C)[1] <- 'Type'
C1 <- gather(C,cell,expr,2:ncol(C))
ggboxplot(C1, x= 'cell', y='expr',
          ylab = "cell", color = 'Type', palette = "jco", merge = "flip", add="jitter")+
  stat_compare_means(size = 4,aes(group= Type),label = 'p.signif',label.x = 1.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5))####  横坐标倾斜45度


##################补充 ssGSEA 推算  ####
setwd("E:/zy/0项目1-40/82.YQ216-10/赵越-2022-01-08-YQ216-10补充/")
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
gene_set<- read.table('../../32.YQ071/08.两种免疫推算/ssGSEA/mmc3.txt',sep = '\t',header = T)
##读取已经下载好的免疫细胞和对应基因列表，来源见文献附件
gene_set<-gene_set[, 1:2]#选取特异基因和对应的免疫细胞两行
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
exprSet=read.table("../00.data/GSE65391/GSE65391_series_matrix.id.uniq.group.sort.txt",
                   header=TRUE,row.names=1,check.names = FALSE)
gsva_matrix <- gsva(as.matrix(exprSet), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
write.table(gsva_matrix,file = 'ssGSEA推算.csv',quote = F,sep = ',')

##########   差异热图
library(RColorBrewer)
library(pheatmap)
gene.a <- read.table("heatmap.txt",header=T,sep="\t",check.names = FALSE,row.names = 1)
annotation_col=read.table("../01.DEG/GSE65391/anno.txt",header=TRUE,
                          row.names=1,check.names = FALSE,sep = '\t')
#gene.a <- log2(gene.a+1)
# gene.b <- t(scale(t(gene.a)))
# gene.b[which(gene.b > quantile(gene.b,0.99))] <- quantile(gene.b,0.99)
# gene.b[which(gene.b < quantile(gene.b,0.01))] <- quantile(gene.b,0.01)
# rownames(gene.b) <- gene.a$id
pheatmap(gsva_matrix,clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         fontsize = 10,clustering_method = "ward.D",
         annotation_col =annotation_col,
         show_rownames = T, #不显示行名
         show_colnames = F,
         treeheight_row = 50,treeheight_col = 50,
         color = colorRampPalette(rev(brewer.pal(n = 9, name ="RdYlBu")))(100),
         border_color = NA, cluster_cols = FALSE)


############   箱线图
library(dplyr)
library(tidyr)
library(ggpubr)
library(cowplot)
boxdat33 <- read.table('box.txt',sep = '\t',header = T,check.names = F)
colnames(boxdat33)[1] <- 'Type'
boxdat44 <- gather(boxdat33,gene,expr,2:ncol(boxdat33))
pdf("ssGSEA.estimate箱线图.pdf",height=8,width=14)
ggboxplot(boxdat44, x= 'gene', y='expr',
          ylab = "proportion", xlab = "cell", color = 'Type', palette = "jco", merge = "flip",
          add="jitter",add.params = list(size=0.5))+
  stat_compare_means(size = 4,aes(group= Type), label = 'p.signif')+
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))
dev.off()







##################    12.相关性   ####
setwd("E:/JZH/12.相关性/")
library(psych)
expr <- read.table('../12.相关性/诊断基因表达量.txt',header = T,row.names = 1,check.names = F,sep = '\t',)
ssgsea <- read.table('18个差异的免疫细胞.txt',header = T,sep = '\t',
                     check.names = F,row.names = 1)
expr <- t(expr)
genelist <- rownames(expr)
expr <- apply(expr,2,function(x){as.numeric(x)})
rownames(expr) <- genelist

#批量计算相关性
gene <- genelist
immuscore <- function(gene){
  y <- as.numeric(expr[gene,])
  colnames <- colnames(ssgsea)
  do.call(rbind,lapply(colnames, function(x){
    dd  <- corr.test(as.numeric(ssgsea[,x]), y , method="spearman",adjust = "fdr")
    data.frame(gene=gene,immune_cells=x,cor=dd$r,p.value=dd$p )
  }))
}


#批量计算genelist跟免疫浸润相关性的结果
data <- do.call(rbind,lapply(genelist,immuscore))
head(data)
#保存到文件
write.csv(data, "ssGSEA.correlation.csv", quote = F, row.names = F)
#增加一列，区分p值的大小,使用两个ifelse实现三分类
data$pstar <- ifelse(data$p.value < 0.05,
                     ifelse(data$p.value < 0.01,"**","*"),
                     "")
data$pstar[1:20]

#开始画图
#相关性用颜色的不同来表示，相关性的大小用颜色的深浅来反映；
#有差异的把*号打印在热图上
#png('Fig29.correlation.png',width = 500,height = 250)
pdf('correlation.pdf',width = 12,height = 8)
ggplot(data, aes(immune_cells, gene)) +
  geom_tile(aes(fill = cor), colour = "black",size=1)+
  scale_fill_gradient2(low = "#2b8cbe",mid = "white",high = "#e41a1c")+
  geom_text(aes(label=pstar),col ="black",size = 7)+
  theme_minimal()+# 不要背景
  theme(axis.title.x=element_blank(),#不要title
        axis.ticks.x=element_blank(),#不要x轴
        axis.title.y=element_blank(),#不要y轴
        axis.text.x = element_text(angle = 45, hjust = 1,size = 11),# 调整x轴文字
        axis.text.y = element_text(size = 11))+#调整y轴文字
  #调整legen
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","Correlation"))
dev.off()
