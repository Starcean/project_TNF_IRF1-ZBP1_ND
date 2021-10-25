#----------------------------------------GSE85162 AD mouse 58 samp-----------------------------------
## 加载R包
library(GEOquery)
## 下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE85162', destdir=".",getGPL = F)
## 获取ExpressionSet对象，包括的表达矩阵和分组信息
gset=gset[[1]]
## 获取分组信息,点开查阅信息
pdata=pData(gset)
pdata
exprSet=exprs(gset)
exprSet
boxplot(exprSet,outline=FALSE, notch=T, las=2)
fix(exprSet)
nrow(exprSet)
ncol(exprSet)

library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)

exprSet = as.data.frame(exprSet)
fix(exprSet)

#判断是否需要log2转化并在需要情况下自动转化
ex <- exprSet
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)

if (LogC) { ex[which(ex <= 0)] <- NaN
exprSet <- log2(ex)
print("log2 transform finished")}else{print("log2 transform not needed")}

#设置工作目录
setwd("D:\\软件\\科研\\百度云备份\\博士阶段\\备份文件\\老婆\\文章修改\\TNF-IRF1-ZBP1-neuron文章修改\\大数据分析相关性")

dt1 <- t(exprSet)
dt1 <- as.data.frame(dt1)
fix(dt1)
Tnf <- "21926_at"
Irf1 <- "16362_at"
Zbp1 <- "58203_at"

n1 <- which(colnames(dt1) == "21926_at")
n2 <- which(colnames(dt1) == "16362_at")
n3 <- which(colnames(dt1) == "58203_at")
Tnf <- dt1[,n1]
Irf1 <- dt1[,n2]
Zbp1 <- dt1[,n3]

plot(Irf1, Zbp1, main = "Irf1_Zbp1_相关性", xlab = "Irf1", ylab = "Zbp1")
cor.test(Irf1,Zbp1,method = "pearson")

#小鼠AD中结果小结：GSE85162 GPL21382 结果 AD mouse IRF1——ZBP1 相关性p=0.0017 


#----------------------------------------GSE151460 AD mouse 122 samp-----------------------------------



# 设置工作空间
setwd('D:\\软件\\科研\\百度云备份\\博士阶段\\备份文件\\老婆\\文章修改\\TNF-IRF1-ZBP1-neuron文章修改\\大数据分析相关性\\GSE151460_RAW\\GSM4578977_JP212.counts.txt')
# 读取该工作空间下的所有文件名
filenames <- dir()
# 通过正则，获取所有txt结尾的文件名
filenames2 <- grep('\\.txt', filenames, value = TRUE)

# 初始化数据框，用于后面的数据合并，因为第一个，第二个样本分隔符sep不同与之后的，所以单独合并
data <- read.table('GSM4578858_JP001.counts.txt', sep = ',',header = T,row.names = 1)
colnames(data) <- c('ID','value')
fix(data)
data2 <- read.table('GSM4578859_JP002.counts.txt', sep = ',',header = T,row.names = 1)
colnames(data2) <- c('ID','value')
fix(data2)
data3 <- merge(data,data2,by = 'ID')
fix(data3)
#-------------------------------------通过循环完成数据合并多个txt文件为一个csv----------------------------------------------

for (i in filenames2[3:122]){
  # 构造数据路径
  dt <- read.table(file = i, sep = '\t',header = T,row.names = 1)
  dt$ID <- row.names(dt)
  data3 <- merge(data3,dt,by = 'ID')
}
write.csv(data3,'data_merg.csv')
#-----------------------------------------------------------------------------------

data_merg <- read.table('data_merg.csv', sep = ',',header = T,row.names = 1)
fix(data_merg)
data_merg=normalizeBetweenArrays(data_merg)
write.csv(data_merg,'data_merg_normalized.csv')

df <- read.csv('data_merg_normalized.csv', sep = ',',header = T,row.names = 1)
fix(df)
t_df <- t(df)
fix(t_df)
class(t_df)
t_df <- as.data.frame(t_df)
plot(t_df$ENSMUSG00000018899,t_df$ENSMUSG00000027514)
cor.test(t_df$ENSMUSG00000018899,t_df$ENSMUSG00000027514)
library(ggplot2)

sp <- ggplot(t_df,aes(x=ENSMUSG00000018899,y=ENSMUSG00000027514))
sp+geom_point()+stat_smooth(method = lm)
#小鼠AD模型中结果小结-IRF1与ZBP1正相关p=0.047，其他无相关


#RNAseq基因ID转换 Mouse ENSMBL_ID to SYMBOL
library(clusterProfiler)
library('org.Mm.eg.db')
eg = bitr(row.names(data_merg), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Mm.eg.db")
class(eg)
write.csv(eg,'ENSEMBL_SYMBOL.CSV')

#----------------------------------------GSE135057 HD mouse 16 samp-----------------------------------

setwd('D:\\软件\\科研\\百度云备份\\博士阶段\\备份文件\\老婆\\文章修改\\TNF-IRF1-ZBP1-neuron文章修改\\大数据分析相关性\\GSE135057_Processed_mu_HD_16sp')
df <- read.table('GSE135057_Processed.txt',header = T,row.names = 1,sep = '\t')
fix(df)
class(df)
t_df <- t(df)
t_df <- as.data.frame(t_df)
fix(t_df)
Tnf <- "ENSMUSG00000024401"
Irf1 <- "ENSMUSG00000018899"
Zbp1 <- "ENSMUSG00000027514"
t_df$Tnf <- t_df$"ENSMUSG00000024401"
t_df$Irf1 <- t_df$"ENSMUSG00000018899"
t_df$Zbp1 <- t_df$"ENSMUSG00000027514"
dt <- t_df

library(ggplot2)
mytheme <- theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#fbfafa"),
        panel.border = element_rect(colour = "black"),
        axis.title.x = element_text(colour = "black",size = 14),
        axis.text.x = element_text(colour = "black",size = 13),
        axis.title.y = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 13))
#相关性
p <- ggplot(data=dt, aes(x=Tnf, y=Irf1))+geom_point(colour="blue",alpha=.8,size=4)+stat_smooth(method="lm",se=F, colour = "red" )
p+mytheme+ggtitle("GSE135057_HD mouse model")+theme(plot.title = element_text(hjust = 0.5))

cor.test(dt$Irf1,dt$Tnf,method = "spearman")

#小鼠HD模型中结果小结-IRF1与ZBP1正相关p=0.047，其他无相关

#----------------------------------------GSE31458 PD mouse 29 samp-----------------------------------

setwd('D:\\软件\\科研\\百度云备份\\博士阶段\\备份文件\\老婆\\文章修改\\TNF-IRF1-ZBP1-neuron文章修改\\大数据分析相关性\\GSE31458_series_matrix_mu_PD_29sp')
df <- read.table('GSE31458_series_matrix.txt',header = T,row.names = 1,sep = '\t')

library(GEOquery)
## 下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE31458', destdir=".",getGPL = F)
## 获取分组信息,点开查阅信息
pdata=pData(gset)
pdata$ID <- row.names(pdata)

library(limma)
df <- normalizeBetweenArrays(df)
boxplot(df,outline=FALSE, notch=T, las=2)
class(df)
t_df <- t(df)
t_df <- as.data.frame(t_df)

Tnf <- "1419607_at"
Irf1 <- "1448436_a_at"
Zbp1 <- "1419604_at"

dt <- t_df
dt$Tnf <- dt$"1419607_at"
dt$Irf1 <- dt$"1448436_a_at"
dt$Zbp1 <- dt$"1419604_at"

dt$ID <- row.names(dt)
dt1 <- merge(dt,pdata,by = "ID")
dt1 <- within(dt1,{
  group <- NA
  group[grepl("control",dt1$title)] <- "WT" #grep()返回value或索引值，grepl()返回逻辑值
  group[grepl("MPTP",dt1$title)] <- "PD"
})
dt1$group <- factor(dt1$group,levels = c("WT","PD"))
#相关性
library(ggplot2)
p <- ggplot(data=dt, aes(x=Tnf, y=Irf1))+geom_point(colour="blue",alpha=.8,size=4)+stat_smooth(method="lm",se=F, colour = "red" )
p+mytheme+ggtitle("GSE31458_PD mouse model")+theme(plot.title = element_text(hjust = 0.5))
cor.test(dt$Irf1,dt$Zbp1,method = "pearson")#相关性检验，p及r value

#ggplot2绘制箱线图
library(ggplot2)
p <- ggplot(dt1,aes(x=factor(group),y=Zbp1))+geom_violin(trim = F)+geom_boxplot(width=.2,fill=c("#00ba38", "#f8786f"))
p+mytheme+ggtitle("GSE74441_AD mouse model")+theme(plot.title = element_text(hjust = 0.5))

#ZBP1在AD高表达
x <- dt1$Zbp1[which(dt3$group=="WT")]
y <- dt1$Zbp1[which(dt3$group=="APPPS1")]
t.test(x,y)

#小鼠PD模型中结果小结TNF-IRF1-ZBP1两两显著正相关,但三者在PD中均无高表达


#----------------------------------------GSE114517 PD Human 75 samp-----------------------------------
# 设置工作空间
setwd('D:\\软件\\科研\\百度云备份\\博士阶段\\备份文件\\老婆\\文章修改\\TNF-IRF1-ZBP1-neuron文章修改\\大数据分析相关性\\GSE114517_RAW')
# 读取该工作空间下的所有文件名
filenames <- dir()
# 通过正则，获取所有txt结尾的文件名
filenames2 <- grep('\\.txt', filenames, value = TRUE)

# 初始化数据框，用于后面的数据合并，因为第一个，第二个样本分隔符sep不同与之后的，所以单独合并
data <- read.table('GSM3143295_SN1.tabular.txt', sep = '\t',header = T,row.names = 1)
data$ID <- row.names(data)
fix(data)

#-------------------------------------通过循环完成数据合并多个txt文件为一个csv----------------------------------------------

for (i in filenames2[2:75]){
  # 构造数据路径
  dt <- read.table(file = i, sep = '\t',header = T,row.names = 1)
  dt$ID <- row.names(dt)
  data <- merge(data,dt,by = 'ID')
}
write.csv(data,'data_merg.csv')

library(limma) 
data_merg <- read.table('data_merg.csv', sep = ',',header = T,row.names = 1)
fix(data_merg)
ncol(data_merg)
nrow(data_merg)
data_merg=normalizeBetweenArrays(data_merg)
write.csv(data_merg,'data_merg_normalized.csv')

df <- read.csv('data_merg_normalized.csv', sep = ',',header = T,row.names = 1)
fix(df)

t_df <- t(df)
fix(t_df)
class(t_df)
t_df <- as.data.frame(t_df)


plot(t_df$ENSG00000125347,t_df$ENSG00000124256)
cor.test(t_df$ENSG00000125347,t_df$ENSG00000124256)
library(ggplot2)

sp <- ggplot(t_df,aes(x=ENSMUSG00000018899,y=ENSMUSG00000027514))
sp+geom_point()+stat_smooth(method = lm)
#小鼠AD模型中结果小结-IRF1与ZBP1正相关p=0.047，其他无相关


#RNAseq基因ID转换 Mouse ENSMBL_ID to SYMBOL
library(clusterProfiler)
library('org.Hs.eg.db')
eg = bitr(row.names(df), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
class(eg)
write.csv(eg,'ENSEMBL_SYMBOL.CSV')

#-------------------------------------------GSE49036 human PD 28samp------------------------------------------------------

## 加载R包
setwd('D:\\软件\\科研\\百度云备份\\博士阶段\\备份文件\\老婆\\文章修改\\TNF-IRF1-ZBP1-neuron文章修改\\大数据分析相关性\\GSE49036_Hs_PD_28sp')
library(GEOquery)

## 下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE49036', destdir=".",getGPL = F)
## 获取ExpressionSet对象，包括的表达矩阵和分组信息
gset=gset[[1]]
## 获取分组信息,点开查阅信息
pdata=pData(gset)
pdata$ID <- row.names(pdata)
fix(pdata)

exprSet=exprs(gset)
boxplot(exprSet,outline=FALSE, notch=T, las=2)

library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)

exprSet = as.data.frame(exprSet)
fix(exprSet)
write.csv(exprSet,"normalized_df.csv")

dt1 <- t(exprSet)
dt1 <- as.data.frame(dt1)
dt1$ID <- row.names(dt1)
dt2 <- merge(dt1,pdata,by = "ID")

dt2$TNF <- dt2$"207113_s_at"
dt2$IRF1 <- dt2$"238725_at"
dt2$ZBP1 <- dt2$"208087_s_at"


dt2 <- within(dt2,{
  group <- NA
  group[grepl("control",dt2$"disease state:ch1")] <- "NC"
  group[grepl("disease",dt2$"disease state:ch1")] <- "PD"
})
dt2$group

library(ggplot2)
mytheme <- theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#fbfafa"),
        panel.border = element_rect(colour = "black"),
        axis.title.x = element_text(colour = "black",size = 14),
        axis.text.x = element_text(colour = "black",size = 13),
        axis.title.y = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 13))

#IRF1 zbp1相关性
p <- ggplot(data=dt2, aes(x=TNF, y=IRF1))+geom_point(colour="blue",alpha=.8,size=4)+stat_smooth(method="lm",se=F, colour = "red" )
p+mytheme+ggtitle("GSE49036_PD Human")+theme(plot.title = element_text(hjust = 0.5))
cor.test(dt2$ZBP1,dt2$IRF1,method = "spearman")#相关性检验，p及r value

#ggplot2绘制箱线图
library(ggplot2)
p <- ggplot(dt2,aes(x=factor(group),y=ZBP1))+geom_violin(trim = F)+geom_boxplot(width=.2,fill=c("#00ba38", "#f8786f"))
p+mytheme+ggtitle("GSE49036_PD Human")+theme(plot.title = element_text(hjust = 0.5))

#ZBP1在AD高表达
x <- dt2$IRF1[which(dt2$group=="NC")]
y <- dt2$IRF1[which(dt2$group=="PD")]
t.test(x,y,var.equal = F)
#人PD模型中结果小结IRF1-ZBP1两两显著正相关
#-------------------------------------------GSE132903 Human AD 195sp-------------------------------------------------------------
## 加载R包
setwd('D:\\软件\\科研\\百度云备份\\博士阶段\\备份文件\\老婆\\文章修改\\TNF-IRF1-ZBP1-neuron文章修改\\大数据分析相关性\\GSE132903_Hs_AD_195sp')
library(GEOquery)

## 下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE132903', destdir=".",getGPL = F)
## 获取ExpressionSet对象，包括的表达矩阵和分组信息
gset=gset[[1]]
## 获取分组信息,点开查阅信息
pdata=pData(gset)
pdata$ID <- row.names(pdata)
write.csv(pdata,"pdata.csv")
class(pdata)
fix(pdata)
pdata$group <- pdata$"diagnosis:ch1"
table(pdata$group)

#ZBP1 10分位、90分位数值进行分组，GSEA
anno <- read.table("GPL10558_annot.txt",sep = "\t",header = T)
dat <- exprSet
dat <- t(dat)
dat <- as.data.frame(dat)
fix(dat)

low <- quantile(dat$ZBP1,0.9)
dat$ID <- row.names(dat)
dat1 <- merge(dat,pdata,by = "ID")
fix(dat1)
dat1$group <- dat1$"diagnosis:ch1"
dat2 <- dat1[which(dat1$group=="AD"),]
dat2$ZBP1 <- dat2$"ILMN_1765994"

quantile(dat2$ZBP1,0.1)

dat3 <- dat2[which(dat2$ZBP1 <= quantile(dat2$ZBP1,0.1)),]
dat3$level <- rep("low",10)
dat4 <- dat2[which(dat2$ZBP1 >= quantile(dat2$ZBP1,0.9)),]
dat4$level <- rep("high",10)
dat5 <- rbind(dat3,dat4)
write.csv(dat5,"dat5.csv")
fix(dat5)
row.names(dat5) <- dat5$ID
which(colnames(dat5) == "ZBP1")
dat6 <- dat5[1:42181]
dat6$level <- dat5$level
dat6$ID <- NULL
fix(dat6)
dat6$ZBP1 <- NULL
dat7 <- t(dat6)
dat7 <- dat7[1:42179,]
#差异表达
as.data.frame(lapply(dat7, as.numeric))
head(dat7)
boxplot(dat7,outline=FALSE)


dat$ID <- row.names(dat)
dat1 <- subset(dat,ID %in% anno$ID)
dat2 <- subset(anno, ID %in% dat$ID)




exprSet=exprs(gset)
exprSet
boxplot(exprSet,outline=FALSE, notch=T, las=2)
fix(exprSet)
nrow(exprSet)
ncol(exprSet)

library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)

exprSet = as.data.frame(exprSet)
fix(exprSet)
write.csv(exprSet,"normalized_df.csv")
head(exprSet)

dt1 <- t(exprSet)
dt1 <- as.data.frame(dt1)

dt1$TNF <- dt1$"ILMN_1728106"
dt1$IRF1 <- dt1$"ILMN_1708375"
dt1$ZBP1 <- dt1$"ILMN_1765994"

dt1$ID <- row.names(dt1)
dt2 <- merge(dt1,pdata,by = "ID")
table(dt2$"diagnosis:ch1")
dt3 <- dt2[which(dt2$group=="AD"),]
nrow(dt3)
dt3$group
table(dt3$group)
dt3$ZBP1
dt4 <- t(dt3)
fix(dt4)
quantile(dt3$ZBP1,c(0.1,0.9))





library(ggplot2)
mytheme <- theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#fbfafa"),
        panel.border = element_rect(colour = "black"),
        axis.title.x = element_text(colour = "black",size = 14),
        axis.text.x = element_text(colour = "black",size = 13),
        axis.title.y = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 13))

#IRF1 zbp1相关性
p <- ggplot(data=dt2, aes(x=IRF1, y=ZBP1))+geom_point(colour="blue",alpha=.8,size=4)+stat_smooth(method="lm",se=F, colour = "red" )
p+mytheme+ggtitle("GSE132903_AD Human")+theme(plot.title = element_text(hjust = 0.5))
cor.test(dt2$ZBP1,dt2$IRF1,method = "pearson")#相关性检验，p及r value

#表达水平分析
dt2$group <- dt2$"diagnosis:ch1"
dt2$group <- factor(dt2$group,levels = c("ND","AD"))

#ggplot2绘制箱线图
library(ggplot2)
p <- ggplot(dt2,aes(x=factor(group),y=ZBP1))+geom_violin(trim = F)+geom_boxplot(width=.2,fill=c("#00ba38", "#f8786f"))
p+mytheme+ggtitle("GSE132903_AD Human")+theme(plot.title = element_text(hjust = 0.5))

#ZBP1在AD高表达
x <- dt2$ZBP1[which(dt2$group=="ND")]
y <- dt2$ZBP1[which(dt2$group=="AD")]
t.test(x,y)

#人AD模型中结果小结IRF1-ZBP1两两显著正相关，IRF1，ZBP1高表达

#--------------------------------------------------------------------------------------------------

#-----------------------------------------GSE85162 mouse 58 AD sample---------------------------------------------------------

## 加载R包
setwd('D:\\软件\\科研\\百度云备份\\博士阶段\\备份文件\\老婆\\文章修改\\TNF-IRF1-ZBP1-neuron文章修改\\大数据分析相关性\\GSE85162_Mu_AD_58sp')
library(GEOquery)

## 下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE85162', destdir=".",getGPL = F)
## 获取ExpressionSet对象，包括的表达矩阵和分组信息
gset=gset[[1]]
## 获取分组信息,点开查阅信息
pdata=pData(gset)
pdata$ID <- row.names(pdata)
write.csv(pdata,"pdata.csv")
class(pdata)
fix(pdata)

exprSet=exprs(gset)
exprSet
boxplot(exprSet,outline=FALSE, notch=T, las=2)
fix(exprSet)
nrow(exprSet)
ncol(exprSet)

library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)

exprSet = as.data.frame(exprSet)
fix(exprSet)
write.csv(exprSet,"normalized_df.csv")
head(exprSet)

dt1 <- t(exprSet)
dt1 <- as.data.frame(dt1)
dt1[1:5,1:5]
Tnf <- "21926_at"
Irf1 <- "58203_at"
Zbp1 <- "16362_at"


dt1$ID <- row.names(dt1)
dt2 <- merge(dt1,pdata,by = "ID")
dt2$Tnf <- dt2$"21926_at"
dt2$Irf1 <- dt2$"58203_at"
dt2$Zbp1 <- dt2$"16362_at"

ncol(dt1)
ncol(pdata)
ncol(dt2)

plot(dt2$"21926_at", dt2$"16362_at", main = "Tnf_Irf1_相关性", xlab = "Tnf", ylab = "Irf1")
abline(lm(formula = dt2$"21926_at"~dt2$"16362_at",data = dt2))

cor.test(dt2$"21926_at",dt2$"16362_at",method = "kendall")

library(ggplot2)
mytheme <- theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#fbfafa"),
        panel.border = element_rect(colour = "black"),
        axis.title.x = element_text(colour = "black",size = 14),
        axis.text.x = element_text(colour = "black",size = 13),
        axis.title.y = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 13))



#IRF1 zbp1相关性
p <- ggplot(data=dt2, aes(x=Irf1, y=Zbp1))+geom_point(colour="blue",alpha=.8,size=4)+stat_smooth(method="lm",se=F, colour = "red" )
p+mytheme+ggtitle("GSE85162_AD mouse model")+theme(plot.title = element_text(hjust = 0.5))

model <- lm(dt2$"58203_at"~dt2$"16362_at",model = T)
summary(model)

#表达水平分析
x1 <- dt2$"21926_at"
x2 <- dt2$Irf1
x3 <- dt2$Zbp1
f <- dt2$"genotype:ch1"
f <- factor(f,levels = c("WT","TG"))
dt2$group <- factor(f,levels = c("WT","TG"))
boxplot(Zbp1~f,dt2,col=c("#00ba38", "#f8786f"))

#ggplot2绘制箱线图
library(ggplot2)
p <- ggplot(dt2,aes(x=factor(group),y=Irf1))+geom_violin(trim = F)+geom_boxplot(width=.2,fill=c("#00ba38", "#f8786f"))
p+mytheme+ggtitle("GSE85162_AD mouse model")+theme(plot.title = element_text(hjust = 0.5))

#ZBP1在AD高表达
x <- dt2$Irf1[which(dt2$"genotype:ch1"=="WT")]
y <- dt2$Irf1[which(dt2$"genotype:ch1"=="TG")]
t.test(x,y)

#小鼠AD模型中结果小结IRF1-ZBP1两两显著正相关，IRF1，ZBP1高表达,Tnf无改变

#-----------------------------------------GSE74441 mouse 260 AD sample---------------------------------------------------------

## 加载R包
setwd('D:\\软件\\科研\\百度云备份\\博士阶段\\备份文件\\老婆\\文章修改\\TNF-IRF1-ZBP1-neuron文章修改\\大数据分析相关性\\GSE74441_Mu_AD_260')
library(GEOquery)

## 下载数据，如果文件夹中有会直接读入
gset = getGEO('GSE74441', destdir=".",getGPL = F)
## 获取ExpressionSet对象，包括的表达矩阵和分组信息
gset=gset[[1]]
## 获取分组信息,点开查阅信息
pdata=pData(gset)
lev <- pdata$"treatment:ch1"
table(factor(lev))

pdata$ID <- row.names(pdata)
write.csv(pdata,"pdata.csv")

pdata_val <- pdata[-which(pdata$"treatment:ch1"=="CD33 ASO"),] #去除CD33 knock-down组
pdata_val$ID <- row.names(pdata_val)#新增ID变量
fix(pdata_val)

exprSet=exprs(gset)
boxplot(exprSet,outline=FALSE, notch=T, las=2)

library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)
exprSet = as.data.frame(exprSet)

dt1 <- t(exprSet) #转置，行名标本编号，变量名探针编号
dt1 <- as.data.frame(dt1)
dt1$ID <- row.names(dt1)#新增ID变量

#探针注释
Tnf <- "ILMN_2899863"
Irf1 <- "ILMN_2649068"
Zbp1 <- "ILMN_2879614"

dt2 <- subset(dt1,ID %in% pdata_val$ID)#筛选出不含CD33 knock-down的样本表达谱
dt3 <- merge(dt2,pdata_val,by = "ID")
dt3$Tnf <- dt3$"ILMN_2899863"
dt3$Irf1 <- dt3$"ILMN_2649068"
dt3$Zbp1 <- dt3$"ILMN_2879614"

dt4 <- dt3[which(dt3$'treatment:ch1'=="PBS"),]

library(ggplot2)
mytheme <- theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#fbfafa"),
        panel.border = element_rect(colour = "black"),
        axis.title.x = element_text(colour = "black",size = 14),
        axis.text.x = element_text(colour = "black",size = 13),
        axis.title.y = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 13))

#IRF1 zbp1相关性
p <- ggplot(data=dt3, aes(x=Irf1, y=Zbp1))+geom_point(colour="blue",alpha=.8,size=4)+stat_smooth(method="lm",se=F, colour = "red" )
p+mytheme+ggtitle("GSE74441_AD mouse model")+theme(plot.title = element_text(hjust = 0.5))

cor.test(dt3$Zbp1,dt3$Irf1,method = "pearson")#相关性检验，p及r value
model <- lm(dt3$Tnf~dt3$Zbp1)
summary(model)

#表达水平分析
a <- dt3$"genotype:ch1"
b <- dt3$"strain:ch1"
length(a)
length(b)
c <- c(a[1:128],rep("APPPS1",45))
table(c)
dt3$group <- c
dt3$group <- factor(dt3$group,levels = c("WT","APPPS1"))

#ggplot2绘制箱线图
library(ggplot2)
p <- ggplot(dt3,aes(x=factor(group),y=Tnf))+geom_violin(trim = F)+geom_boxplot(width=.2,fill=c("#00ba38", "#f8786f"))
p+mytheme+ggtitle("GSE74441_AD mouse model")+theme(plot.title = element_text(hjust = 0.5))

#ZBP1在AD高表达
x <- dt3$Zbp1[which(dt3$group=="WT")]
y <- dt3$Zbp1[which(dt3$group=="APPPS1")]
t.test(x,y)

#小结：小鼠AD中Tnf1-Irf1-Zbp1显著相关，单Zbp1无显著上调

#------------------------------------------AC infection model----------------------------------------------
setwd('D:\\软件\\科研\\百度云备份\\博士阶段\\备份文件\\老婆\\文章修改\\TNF-IRF1-ZBP1-neuron文章修改\\大数据分析相关性\\AC_mouse_Tnf_inhibitor')
dt <- read.table("TNFZBP1.txt",header = T,sep = "\t")
fix(dt)

library(ggplot2)
mytheme <- theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill="#fbfafa"),
        panel.border = element_rect(colour = "black"),
        axis.title.x = element_text(colour = "black",size = 14),
        axis.text.x = element_text(colour = "black",size = 13),
        axis.title.y = element_text(colour = "black",size = 14),
        axis.text.y = element_text(colour = "black",size = 13))


p <- ggplot(dt,aes(x=TNF,y=ZBP1,colour=Group))+geom_point(size=4)
p+geom_abline(intercept = 0.661,slope = 1.495,size=0.75,colour="red")+mytheme+ggtitle("AC infection mouse model")+theme(plot.title = element_text(hjust = 0.5))

cor.test(dt$TNF,dt$ZBP1,method = "pearson")#相关性检验，p及r value

































#ZBP1 COIP 富集分析






























