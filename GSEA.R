#GSEA分析ZBP1在小鼠及人ND中富集的信号通路
#加载包
library(clusterProfiler)
library(org.Mm.eg.db)

#分析GSE74441
setwd("D:\\软件\\科研\\百度云备份\\博士阶段\\备份文件\\老婆\\文章修改\\TNF-IRF1-ZBP1-neuron文章修改\\GSEA\\GSE74441_Mu_AD_260")
d1 <- read.table("GSE74441_series_matrix.txt",sep = "\t",header = T,row.names = 1)
row.names(d1)[1:10]
d1$ID <- row.names(d1)
fix(d1)

#读入整理好的GPL文件
d2 <- read.table("GPL6885.txt",sep = "\t",header = T)
fix(d2)

nrow(d1)
nrow(d2)

#分别取公共ID，子集
d3 <- subset(d1,ID %in% d2$ID)
d4 <- subset(d2,ID %in% d1$ID)

#合并数据框，生成包含ID，SYMBOL，和GENE ID的数据框
d5 <- merge(d3,d4,by = "ID") #58203(Gene_ID of Zbp1)
fix(d5)
#删除ID列
d5$ID <- NULL

#删除GENE ID为空的行
d6 <- d5[-which(is.na(d5$Entrez_Gene_ID)), ]
nrow(d5)
nrow(d6)

#删除GENE ID列，只保留SYMBOL作为标识变量
d6$Entrez_Gene_ID <- NULL

#转置数据框
d7 <- as.data.frame(t(d6)) 

#定义变量名为SYMBOL，注意提取数据框指定行返回的仍为数据框，需要翻转再取列生成变量
colnames(d7) <- t(d7[261,])[,1]

#去除SYMBOL行
d8 <- d7[1:260,]
fix(d8)

#写出d8，变量：，行名：。
write.csv(d8,"d8.csv")

#重新读入d8
d8 <- read.csv("d8.csv",header = T,row.names = 1)

d8 <- within(d8,{
  Zbp1level <- NULL
  Zbp1level[Zbp1>=quantile(d8$Zbp1,probs = 0.9)] <- "High"
  Zbp1level[Zbp1<=quantile(d8$Zbp1,probs = 0.1)] <- "Low"
  Zbp1level[Zbp1<quantile(d8$Zbp1,probs = 0.9)&Zbp1>quantile(d8$Zbp1,probs = 0.1)] <- "Median"
})  

table(d8$Zbp1level) #经计算大于90分位的样本26个，小于10分位的样本也是26个

#提取Zbp1表达水平大于90分位以及小于10分位的样本
d9 <- d8[(d8$Zbp1level=="High" | d8$Zbp1level=="Low"),]
write.csv(d9,"d9.CSV")

#Zbp1高低分组间,差异表达基因计算
library(limma)
library(impute)

#excel修改d9后读入d9，行为基因名，列为分组
rt <- read.csv("d9.csv",header = T)

#删除包含NA的所有行
rt <- rt[!(rowSums(is.na(rt))),]

#删除含Rik的行
#grep()函数按照正则表达式的pattern进行查找对象，并将对象中的将满足正则表达式pattern的元素下标返回来
n <- grep(pattern = 'Rik',x=rt$Zbp1level)
rt <- rt[-n,]

#转为矩阵，且计算缺失值，填补缺失值
rt=as.matrix(rt)
rownames(rt)=rt[,1]
rt=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
exp=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

#impute missing expression data
mat=impute.knn(exp)
rt=mat$data

#标准化处理，去除批次效应
#normalize
pdf(file="rawall.pdf")
boxplot(rt,col = "blue",xaxt = "n",outline = F)
dev.off()
rt=normalizeBetweenArrays(as.matrix(rt))
pdf(file="normalizedall.pdf")
boxplot(rt,col = "red",xaxt = "n",outline = F)
dev.off()

#差异基因计算
class <- c(rep(c("N","T"),each=26))
design <- model.matrix(~0+factor(class))
colnames(design) <- c("N","T")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(T-N,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=10002)
write.table(allDiff,file="all_DEGs1.xls",sep="\t",quote=F)

#读入excel修改好的DEGs
#读入平台文件，并将DEGs中的Symbol转为Gene ID
dt1 <- read.csv("DEGs.csv",header = T)

#Symbol转为Gene ID
library(clusterProfiler)
library(org.Mm.eg.db)
x <- dt1$SYMBOL
dt2 <- bitr(x, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
colnames(dt1) <- c("SYMBOL","FC")

#取子集
dt3 <- subset(dt1,SYMBOL%in%dt2$SYMBOL)
dt4 <- subset(dt2,SYMBOL%in%dt1$SYMBOL)

#合并数据集获取FC及GENE ID 两列数据框
dt5 <- merge(dt3,dt4,by="SYMBOL")
dt5$SYMBOL <- NULL
colnames(dt5) <- c("FC","ID")
fix(dt5)
#获取GENE list用于GSEA分析
#文件为整理好的第一列为FC，第二列为ID,ID为geneID，
geneList = dt5[,1]
names(geneList) = as.character(dt5[,2])
geneList1 = sort(geneList, decreasing = TRUE)
fix(geneList1)
#GSEA对KEGG分析
KK <- gseKEGG(geneList = geneList1,organism = "mmu",nPerm = 100,minGSSize = 5,pvalueCutoff = 1,verbose = FALSE)
write.csv(KK,"KKEGG.csv")
#GSEA对KEGG分析后提取指定基因集
gseaplot(KK,geneSetID = "mmu04217")
gseaplot(KK,geneSetID = "mmu04210")
gseaplot(KK,geneSetID = "mmu04217")

#-----------------------------------------------------------------------------------------------
#分析GSE31458
setwd("D:\\软件\\科研\\百度云备份\\博士阶段\\备份文件\\老婆\\文章修改\\TNF-IRF1-ZBP1-neuron文章修改\\GSEA\\GSE31458_Mu_PD_29")
#读入探针表达矩阵
d1 <- read.table("GSE31458_series_matrix.txt",sep = "\t",header = T,row.names = 1)
d1$ID <- row.names(d1)
#删除最后一行
d1 <- d1[-22691,]

#读入整理过的GPL文件
d2 <- read.table("GPL8321.txt",sep = "\t",header = T,na.strings=" ")

#分别取公共ID，子集
d3 <- subset(d1,ID %in% d2$ID)
d4 <- subset(d2,ID %in% d1$ID)

#合并数据框，生成包含ID，SYMBOL，和GENE ID的数据框
d5 <- merge(d3,d4,by = "ID") #58203(Gene_ID of Zbp1)

#删除ID列
d5$ID <- NULL
fix(d5)

#删除GENE ID为空的行
#删除GENE ID为空的行
d6 <-  d5
nrow(d5)
nrow(d6)

#删除GENE ID列，只保留SYMBOL作为标识变量
d6$ENTREZ_GENE_ID <- NULL

#转置数据框
d7 <- as.data.frame(t(d6)) 
fix(d7)

#定义变量名为SYMBOL，注意提取数据框指定行返回的仍为数据框，需要翻转再取列生成变量
colnames(d7) <- t(d7[30,])[,1]

#去除SYMBOL行
d8 <- d7[1:29,]
fix(d8)

summary(d8)

#-----------------------------------------------------------------------------------------------
#分析GSE49036
setwd("D:\\软件\\科研\\百度云备份\\博士阶段\\备份文件\\老婆\\文章修改\\TNF-IRF1-ZBP1-neuron文章修改\\GSEA\\GSE31458_Mu_PD_29")

#读取基因表达文件
probe_exp<-read.table("GSE31458_series_matrix.txt",header=T,sep="\t",row.names=1)
#读取探针文件
probeid_geneid<-read.table("GPL8321.txt",header=T,sep="\t")
fix(probeid_geneid)
probe_name<-rownames(probe_exp)
#probe进行匹配
loc<-match(probeid_geneid[,1],probe_name)
#确定能匹配上的probe表达值
probe_exp<-probe_exp[loc,]
#每个probeid对应的geneid
raw_geneid<-as.numeric(as.matrix(probeid_geneid[,3]))
#找出有geneid的probeid并建立索引
index<-which(!is.na(raw_geneid))
#提取有geneid的probe
geneid<-raw_geneid[index]
#找到每个geneid的表达值
exp_matrix<-probe_exp[index,]
geneidfactor<-factor(geneid)
#多个探针对应1个基因的情况，取平均值
gene_exp_matrix<-apply(exp_matrix,2,function(x) tapply(x,geneidfactor,mean))
#geneid作为行名
rownames(gene_exp_matrix)<-levels(geneidfactor)
geneid<-rownames(gene_exp_matrix)
gene_exp_matrix2<-cbind(geneid,gene_exp_matrix)
write.table(gene_exp_matrix2,file="geneid.exprs.txt",sep='\t',quote=F,row.names=F)
#将gene id 转换为gene symbol
loc<-match(rownames(gene_exp_matrix),probeid_geneid[,3])
rownames(gene_exp_matrix)=probeid_geneid[loc,2]
genesymbol<-rownames(gene_exp_matrix)
gene_exp_matrix3<-cbind(genesymbol,gene_exp_matrix)
write.table(gene_exp_matrix3,file="genesyb.exprs.txt",sep='\t',quote=F,row.names=F)


#转为矩阵，且计算缺失值，填补缺失值
rt=as.matrix(gene_exp_matrix3)

rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
exp=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)

library(limma)
library(impute)
#impute missing expression data
mat=impute.knn(exp)
rt=mat$data

#标准化处理，去除批次效应
#normalize
pdf(file="rawall.pdf")
boxplot(rt,col = "blue",xaxt = "n",outline = F)
dev.off()
rt=normalizeBetweenArrays(as.matrix(rt))
pdf(file="normalizedall.pdf")
boxplot(rt,col = "red",xaxt = "n",outline = F)
dev.off()



#
dat <- read.table("genesyb.exprs.txt",header=T,sep="\t")
dat <- as.data.frame(t(dat))
fix(dat)
colnames(dat) <- t(dat[1,])[,1]
dat <- dat[-1,]
write.csv(dat,"genesyb_exp.csv")

rt <- read.csv("genesyb_exp.csv",header = T)
name <- rt[,1] 
rt <- rt[,-1]
rt <- t(rt)
colnames(rt) <- name
fix(rt)
rt[which(row.names(rt)=="Zbp1"),]

class <- c(rep(c("high","low"),each=7))
design <- model.matrix(~0+factor(class))
colnames(design) <- c("high","low")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(high-low,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=13046)
write.table(allDiff,file="all_DEGs1.xls",sep="\t",quote=F)
nrow(rt)


#Symbol转为Gene ID
library(clusterProfiler)
library(org.Mm.eg.db)
dt1 <- read.csv("forGSEA.csv",header = T)
fix(dt1)
x <- dt1$SYMBOL

dt2 <- bitr(x, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Mm.eg.db")
class(dt2)
colnames(dt2)
#取子集
dt3 <- subset(dt1,SYMBOL%in%dt2$SYMBOL)
dt4 <- subset(dt2,SYMBOL%in%dt1$SYMBOL)

#合并数据集获取FC及GENE ID 两列数据框
dt5 <- merge(dt3,dt4,by="SYMBOL")
dt5$SYMBOL <- NULL
colnames(dt5) <- c("FC","ID")
fix(dt5)
#获取GENE list用于GSEA分析
#文件为整理好的第一列为FC，第二列为ID,ID为geneID，
geneList = dt5[,1]
names(geneList) = as.character(dt5[,2])
geneList1 = sort(geneList, decreasing = TRUE)
fix(geneList1)
#GSEA对KEGG分析
KK <- gseKEGG(geneList = geneList1,organism = "mmu",nPerm = 500,minGSSize = 20,pvalueCutoff = 1,verbose = FALSE)
write.csv(KK,"KKEGG.csv")
#GSEA对KEGG分析后提取指定基因集
gseaplot(KK,geneSetID = "mmu04210")
gseaplot(KK,geneSetID = "mmu04140")
gseaplot(KK,geneSetID = "mmu04145")























