#-----------------------------------------GSE74441 mouse 260 AD sample---------------------------------------------------------

## set workspace
setwd('location')
library(GEOquery)

## download data set
gset = getGEO('GSE74441', destdir=".",getGPL = F)
## obtain expression matrix and group information
gset=gset[[1]]
pdata=pData(gset)
lev <- pdata$"treatment:ch1"
pdata$ID <- row.names(pdata)
write.csv(pdata,"pdata.csv")

pdata_val <- pdata[-which(pdata$"treatment:ch1"=="CD33 ASO"),] #去除CD33 knock-down组
pdata_val$ID <- row.names(pdata_val)#新增ID变量

exprSet=exprs(gset)
boxplot(exprSet,outline=FALSE, notch=T, las=2)

#normalize
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)
exprSet = as.data.frame(exprSet)

dt1 <- t(exprSet) 
dt1 <- as.data.frame(dt1)
dt1$ID <- row.names(dt1)

#probe annotation
Tnf <- "ILMN_2899863"
Irf1 <- "ILMN_2649068"
Zbp1 <- "ILMN_2879614"

#subset needed data
dt2 <- subset(dt1,ID %in% pdata_val$ID)
dt3 <- merge(dt2,pdata_val,by = "ID")
dt3$Tnf <- dt3$"ILMN_2899863"
dt3$Irf1 <- dt3$"ILMN_2649068"
dt3$Zbp1 <- dt3$"ILMN_2879614"

dt4 <- dt3[which(dt3$'treatment:ch1'=="PBS"),]

#set plot theme
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

#Correlation between IRF1 and zbp1
p <- ggplot(data=dt3, aes(x=Irf1, y=Zbp1))+geom_point(colour="blue",alpha=.8,size=4)+stat_smooth(method="lm",se=F, colour = "red" )
p+mytheme+ggtitle("GSE74441_AD mouse model")+theme(plot.title = element_text(hjust = 0.5))

cor.test(dt3$Zbp1,dt3$Irf1,method = "pearson")
model <- lm(dt3$Tnf~dt3$Zbp1)
summary(model)

#expression level
a <- dt3$"genotype:ch1"
b <- dt3$"strain:ch1"
c <- c(a[1:128],rep("APPPS1",45))
dt3$group <- c
dt3$group <- factor(dt3$group,levels = c("WT","APPPS1"))

#plot boxplot
library(ggplot2)
p <- ggplot(dt3,aes(x=factor(group),y=Tnf))+geom_violin(trim = F)+geom_boxplot(width=.2,fill=c("#00ba38", "#f8786f"))
p+mytheme+ggtitle("GSE74441_AD mouse model")+theme(plot.title = element_text(hjust = 0.5))

#ZBP1 expression
x <- dt3$Zbp1[which(dt3$group=="WT")]
y <- dt3$Zbp1[which(dt3$group=="APPPS1")]
t.test(x,y)

#----------------------------------------GSE31458 PD mouse 29 samp-----------------------------------
#Set workspace
setwd('location')
#Reading file
df <- read.table('GSE31458_series_matrix.txt',header = T,row.names = 1,sep = '\t')

#loading package
library(GEOquery)
## Downloading data set
gset = getGEO('GSE31458', destdir=".",getGPL = F)
## Obatain group information
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
  group[grepl("control",dt1$title)] <- "WT" 
  group[grepl("MPTP",dt1$title)] <- "PD"
})
dt1$group <- factor(dt1$group,levels = c("WT","PD"))

#correlation calculate
library(ggplot2)
p <- ggplot(data=dt, aes(x=Tnf, y=Irf1))+geom_point(colour="blue",alpha=.8,size=4)+stat_smooth(method="lm",se=F, colour = "red" )
p+mytheme+ggtitle("GSE31458_PD mouse model")+theme(plot.title = element_text(hjust = 0.5))
cor.test(dt$Irf1,dt$Zbp1,method = "pearson")

#plot
library(ggplot2)
p <- ggplot(dt1,aes(x=factor(group),y=Zbp1))+geom_violin(trim = F)+geom_boxplot(width=.2,fill=c("#00ba38", "#f8786f"))
p+mytheme+ggtitle("GSE74441_AD mouse model")+theme(plot.title = element_text(hjust = 0.5))

#ZBP1 expression
x <- dt1$Zbp1[which(dt3$group=="WT")]
y <- dt1$Zbp1[which(dt3$group=="APPPS1")]
t.test(x,y)

#-------------------------------------------GSE49036 human PD 28samp------------------------------------------------------

## set workspace
setwd('location')

#load package
library(GEOquery)

## Download data set
gset = getGEO('GSE49036', destdir=".",getGPL = F)
## Obtain expression of matrix and group information
gset=gset[[1]]
pdata=pData(gset)
pdata$ID <- row.names(pdata)

exprSet=exprs(gset)
boxplot(exprSet,outline=FALSE, notch=T, las=2)

library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)
exprSet = as.data.frame(exprSet)
write.csv(exprSet,"normalized_df.csv")

dt1 <- t(exprSet)
dt1 <- as.data.frame(dt1)
dt1$ID <- row.names(dt1)
dt2 <- merge(dt1,pdata,by = "ID")

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

#correaltion between IRF1 and zbp1
p <- ggplot(data=dt2, aes(x=IRF1, y=ZBP1))+geom_point(colour="blue",alpha=.8,size=4)+stat_smooth(method="lm",se=F, colour = "red" )
p+mytheme+ggtitle("GSE49036_PD Human")+theme(plot.title = element_text(hjust = 0.5))
cor.test(dt2$ZBP1,dt2$IRF1,method = "spearman")

#plot
library(ggplot2)
p <- ggplot(dt2,aes(x=factor(group),y=ZBP1))+geom_violin(trim = F)+geom_boxplot(width=.2,fill=c("#00ba38", "#f8786f"))
p+mytheme+ggtitle("GSE49036_PD Human")+theme(plot.title = element_text(hjust = 0.5))

#ZBP1expression
x <- dt2$ZBP1[which(dt2$group=="NC")]
y <- dt2$ZBP1[which(dt2$group=="PD")]
t.test(x,y,var.equal = F)

#-------------------------------------------GSE132903 Human AD 195sp-------------------------------------------------------------
## set workspace
setwd('location')

#load package
library(GEOquery)

## download data set
gset = getGEO('GSE132903', destdir=".",getGPL = F)
## obtain expression of matrix and group information
gset=gset[[1]]
pdata=pData(gset)
pdata$ID <- row.names(pdata)
write.csv(pdata,"pdata.csv")
pdata$group <- pdata$"diagnosis:ch1"
table(pdata$group)

#subsetting data forGSEA
anno <- read.table("GPL10558_annot.txt",sep = "\t",header = T)
dat <- exprSet
dat <- t(dat)
dat <- as.data.frame(dat)
low <- quantile(dat$ZBP1,0.9)
dat$ID <- row.names(dat)
dat1 <- merge(dat,pdata,by = "ID")
dat1$group <- dat1$"diagnosis:ch1"
dat2 <- dat1[which(dat1$group=="AD"),]
dat2$ZBP1 <- dat2$"ILMN_1765994"
dat3 <- dat2[which(dat2$ZBP1 <= quantile(dat2$ZBP1,0.1)),]
dat3$level <- rep("low",10)
dat4 <- dat2[which(dat2$ZBP1 >= quantile(dat2$ZBP1,0.9)),]
dat4$level <- rep("high",10)
dat5 <- rbind(dat3,dat4)
write.csv(dat5,"dat5.csv")
row.names(dat5) <- dat5$ID
which(colnames(dat5) == "ZBP1")
dat6 <- dat5[1:42181]
dat6$level <- dat5$level
dat6$ID <- NULL
dat6$ZBP1 <- NULL
dat7 <- t(dat6)
dat7 <- dat7[1:42179,]

#normalize
exprSet=exprs(gset)
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T, las=2)

exprSet = as.data.frame(exprSet)
write.csv(exprSet,"normalized_df.csv")

dt1 <- t(exprSet)
dt1 <- as.data.frame(dt1)

dt1$TNF <- dt1$"ILMN_1728106"
dt1$IRF1 <- dt1$"ILMN_1708375"
dt1$ZBP1 <- dt1$"ILMN_1765994"

dt1$ID <- row.names(dt1)
dt2 <- merge(dt1,pdata,by = "ID")
dt3 <- dt2[which(dt2$group=="AD"),]
dt4 <- t(dt3)

#set plot theme
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

#correaltion between IRF1 and zbp1
p <- ggplot(data=dt2, aes(x=IRF1, y=ZBP1))+geom_point(colour="blue",alpha=.8,size=4)+stat_smooth(method="lm",se=F, colour = "red" )
p+mytheme+ggtitle("GSE132903_AD Human")+theme(plot.title = element_text(hjust = 0.5))
cor.test(dt2$ZBP1,dt2$IRF1,method = "pearson")

#expression
dt2$group <- dt2$"diagnosis:ch1"
dt2$group <- factor(dt2$group,levels = c("ND","AD"))

#plot
library(ggplot2)
p <- ggplot(dt2,aes(x=factor(group),y=ZBP1))+geom_violin(trim = F)+geom_boxplot(width=.2,fill=c("#00ba38", "#f8786f"))
p+mytheme+ggtitle("GSE132903_AD Human")+theme(plot.title = element_text(hjust = 0.5))

#ZBP1 expression
x <- dt2$ZBP1[which(dt2$group=="ND")]
y <- dt2$ZBP1[which(dt2$group=="AD")]
t.test(x,y)

#----------------------------------------GSE135057 HD mouse 16 samp-----------------------------------
#Set workspace
setwd('location')
#Reading file
df <- read.table('GSE135057_Processed.txt',header = T,row.names = 1,sep = '\t')
t_df <- t(df)
t_df <- as.data.frame(t_df)
Irf1 <- "ENSMUSG00000018899"
Zbp1 <- "ENSMUSG00000027514"
t_df$Irf1 <- t_df$"ENSMUSG00000018899"
t_df$Zbp1 <- t_df$"ENSMUSG00000027514"
dt <- t_df

#Plot
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

#Correlation calculating
p <- ggplot(data=dt, aes(x=Irf1, y=Zbp1))+geom_point(colour="blue",alpha=.8,size=4)+stat_smooth(method="lm",se=F, colour = "red" )
p+mytheme+ggtitle("GSE135057_HD mouse model")+theme(plot.title = element_text(hjust = 0.5))
cor.test(dt$Irf1,dt$Tnf,method = "spearman")

#------------------------------------------AC infection model----------------------------------------------
setwd('location')
dt <- read.table("expression matrix.txt",header = T,sep = "\t")

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

cor.test(dt$TNF,dt$ZBP1,method = "pearson")



