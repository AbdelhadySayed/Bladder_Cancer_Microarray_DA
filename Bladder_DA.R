####load libraries ########  
library(readxl)
library(readr)
library(ggfortify)
library(matrixStats)
library(rgl)
library(scatterplot3d)
library(ComplexHeatmap)
library(genefilter)
library(multtest)
library(Biobase)
library(affy)
library(affyPLM)
library(hgu133a.db)
library(hgu133plus2.db)
library(hgu133acdf)
library(hgu133plus2cdf)


### load the CEL files
celFilesDirectory="GSE7476"
cels = list.files(celFilesDirectory, pattern = "CEL")
affyData <- ReadAffy(celfile.path=celFilesDirectory)
affyData



### Data Exploration
class(affyData)
sampleNames(affyData)
featureNames(affyData)
head(featureNames(affyData))
tail(featureNames(affyData))
annotation(affyData)
dim(affyData)


# see how the RAW expression look like without processing : notice the large values, exprs:to extract the expression matrix
head(exprs(affyData))
head(Biobase::exprs(affyData))
View(Biobase::exprs(affyData))

exp.raw= Biobase::exprs(affyData)


# read metadata
meta <- read_excel("GSE7476/phenotable.xlsx")
meta

meta$sample.type2 = unlist(sapply(strsplit(meta$Sample.Type, "_"),c)[2,])
meta

# PCA
pca_rw <- prcomp(t(exp.raw), scale. = TRUE)
autoplot(pca_rw, data = meta, colour = 'sample.type2')

scores = as.data.frame(pca_rw$x)
head(scores[1:4])

colors <- c("blue", "red", "brown", "green") # Number of color according to the number of groups
colors <- colors[as.numeric(meta$sample.type2)]
meta$sample.type2 = as.factor(meta$sample.type2)

#s3d <- scatterplot3d(exp.pca$loadings[,1:3], pch=".",col.axis="blue", col.grid="lightblue")

s3d <-scatterplot3d(pca_rw$x[, 1], pca_rw$x[, 2],pca_rw$x[, 3],xlab="PC1",ylab="PC2", zlab="PC3", pch = 16,color=colors)
legend("right", legend = levels(meta$sample.type2),
  col =  c("blue", "red", "brown", "green") , pch = 16)


# another way to explore expressions in the first 3 genes/probes in the first 5 columns
Biobase::exprs(affyData)[1:3 , 1:5]

# 1- histogram, length of sampleNames=12 from above
hist(affyData,main = "Histogram affyData")
cols=seq(1: length( sampleNames(affyData) ))
legend("topright",sampleNames(affyData),col=cols,lty=1,lwd=3,cex=0.3,ncol = 2)####

boxplot(affyData,main = "Box Plot ",col=seq(1: length( sampleNames(affyData) )), las=2)



### Microarray data preprocessing
# data pre-processing in one step : life is so easy !
# threestep (background correction, normalization, summarization from probes to probesets)
eset = threestep(affyData,
                 background.method = "IdealMM",
                 normalize.method = "quantile",
                 summary.method = "average.log")

#eset2 = expresso(affyData, bgcorrect.method="rma", 
                           #normalize.method="constant",pmcorrect.method="pmonly",
                           #summary.method="avgdiff")
# the expresso function doesn't do log transformation.
#don't forget to do it yourrself... .check the ranges 
##View(exprs(eset2))
##range(exprs(eset2))

View(exprs(eset))
range(exprs(eset))

# threestep (background correction, normalization, summarization from probes to probesets)
eset2 = threestep(affyData,
                 background.method = "MAS",
                 normalize.method = "quantile",
                 summary.method = "average.log")

#eset2 = expresso(affyData, bgcorrect.method="rma", 
                           #normalize.method="constant",pmcorrect.method="pmonly",
                           #summary.method="avgdiff")
# the expresso function doesn't do log transformation.
#don't forget to do it yourrself... .check the ranges 
##View(exprs(eset2))
##range(exprs(eset2))

View(exprs(eset2))
range(exprs(eset2))

## view data after preprocessing
hist(eset2,main = "Histogram affyData")
boxplot(eset2,main = "Box Plot GDS596",col=seq(1:23), las=2)
exp=exprs(eset2)
View(exp)


data=cbind(probe_id=rownames(exp),exp)
colnames(data)
head(data)


#### mapping the probe ids into gene symbole

ls("package:hgu133plus2.db")
hgu133plus2()
mapper = hgu133plus2SYMBOL
mapper 
map.df = as.data.frame(mapper)
head(map.df)

# merge the two data frames to have the symbole annotation in the data object
data2=merge(data,map.df,by="probe_id",all.x=T)
head(data2)

# do i need the probe id again?  no, then let's drop it
data2=data2[,-1]
# remove nulls : some probes were not mapped to any gene symbol, to remove the rows has no gene symbol
sum(is.na(data2$symbol))
data2=data2[ ! is.na(data2$symbol), ]
# check duplciation of of gene symbols?  
x=duplicated(data2$symbol)
sum(x)

### yes .. why ? probesets?  solutions : aggregation
exp.data=data2[-dim(data2)[2]]
exp.data=apply(exp.data,2, as.numeric)

####remove duplication
exp.data.agg= aggregate(exp.data, by=list(data2$symbol),FUN=mean)
names(exp.data.agg)
dim(exp.data.agg)
# make genes(group.1) as row names
rownames(exp.data.agg)=exp.data.agg$Group.1
# drop group.1 column
exp.data.agg=exp.data.agg[- 1]
# remove .CEL.gz from col names
names(exp.data.agg)= unlist(sapply(strsplit(names(exp.data.agg), "\\."),c)[1,])


# save the data and metadata files
exp=exp.data.agg
head(exp)

#PCA
pca_preprs <- prcomp(t(exp), scale. = TRUE)
autoplot(pca_preprs, data = meta, colour = 'sample.type2')

scores2 = as.data.frame(pca_preprs$x)
head(scores2[1:4])

#3D PCA
colors <- c("blue", "red", "black", "green") # Number of color according to the number of groups
colors <- colors[as.numeric(meta$sample.type2)]
meta$sample.type2 = as.factor(meta$sample.type2)
#s3d <- scatterplot3d(exp.pca$loadings[,1:3], pch=".",col.axis="blue", col.grid="lightblue")
s3d <-scatterplot3d(pca_preprs$x[, 1], pca_preprs$x[, 2],pca_preprs$x[, 3],xlab="PC1",ylab="PC2", zlab="PC3", pch = 16,color=colors)
legend("right", legend = levels(meta$sample.type2),
  col =  c("blue", "red", "black", "green") , pch = 16)

#Boxplot
boxplot(exp,main = "Box Plot GDS596 Preprocessed",col=seq(1:23), las=2)

#Densityplot
library("tidyr")
stak_dat <- exp %>%                          # Apply pivot_longer function
  pivot_longer(colnames(exp)) %>% 
  as.data.frame()
head(stak_dat)

ggplot(stak_dat, aes(x=value, fill=name)) + geom_density(alpha=0.2)



## Grouping of data based on conditions to make differential analysis
groups=unique(meta$sample.type2)
group1=groups[1]
group2=groups[2]
group3=groups[3]
group4=groups[4]

group1.columns=  meta[ meta$sample.type2 ==group1, ]$Sample.Id
group2.columns=  meta[ meta$sample.type2 ==group2,]$Sample.Id
group3.columns=  meta[ meta$sample.type2 ==group3, ]$Sample.Id
group4.columns=  meta[ meta$sample.type2 ==group4,]$Sample.Id

exp1=exp[ , c(group1.columns, group3.columns)]
head(exp1)

#LogFC
lfc=apply(exp1,1, function(x)  mean(x[group3.columns]) -mean(x[group1.columns]))
res=as.data.frame(lfc)
head(res)

#Pval and adjPval
f=factor( c( rep(1, length(group1.columns)) , rep(2, length(group3.columns)) ))
pval=rowttests(as.matrix(exp1),f)$p.value
correctPvalueandReturnAll<-function(tt.pval,method)
{
  mt=mt.rawp2adjp(tt.pval,proc=method)
  adjp=mt$adjp[order(mt$index),]
  return(adjp[,2])
}
adj.pval=correctPvalueandReturnAll(pval,"BH")
res=cbind(lfc,pval,adj.pval)
head(res)


#####  selection criteria for identifying DEGS
res.degs=res[ abs(lfc) > log2(2) & pval <0.05,]  ##### identify DEGs based on both  LFC and the significance level 
res.degs=res.degs[order(res.degs[,3]) ,]
degs.genes= rownames(res.degs)
exp1.degs=exp1[degs.genes,]
head(exp1.degs)
dim(exp1.degs)

#save degs genes
write.table(degs.genes,file = "DEGs.txt",row.names = F,col.names = F,quote = F)#to get txt file

# T1GH vs control
meta.sub=meta[ meta$sample.type2=="control" | meta$sample.type2=="T1GH" ,]
meta.sub

#PCA for DEGS only
pca_deg <- prcomp(t(exp1.degs), scale. = TRUE)
autoplot(pca_deg, data = meta.sub, colour = 'sample.type2')

#3D PCA
colors <- c("blue", "red") # Number of color according to the number of groups
colors <- colors[as.numeric(meta.sub$sample.type2)]
#meta.sub$sample.type2 = as.factor(meta.sub$sample.type2)

#s3d <- scatterplot3d(exp.pca$loadings[,1:3], pch=".",col.axis="blue", col.grid="lightblue")

s3d <-scatterplot3d(pca_deg$x[, 1], pca_deg$x[, 2],pca_deg$x[, 3],xlab="PC1",ylab="PC2", zlab="PC3", pch = 16,color=colors)
legend("right", legend = levels(meta.sub$sample.type2),
  col =  c("blue", "red") , pch = 16)

# Densityplot
stak_dat3 <- exp1.degs %>%                          # Apply pivot_longer function
  pivot_longer(colnames(exp1.degs)) %>% 
  as.data.frame()
head(stak_dat3)
ggplot(stak_dat3, aes(x=value, fill=name)) + geom_density(alpha=0.2)

# boxplot
boxplot(exp1.degs,main = "Box Plot GDS596 DEGs",col=seq(1:23), las=2)

# Assign UP & DOWN regulated genes
res.degs=as.data.frame(res.degs)
res.degs[,"regulation"]="DOWN"
res.degs[ res.degs$lfc > 0 ,  ]$regulation="UP"
res.degs.up=res.degs[res.degs$regulation=="UP" , ]
res.degs.down=res.degs[res.degs$regulation=="DOWN" , ]
res.degs.down

### Volcanoplot

plot(res[,1], -10*log10(res[,2]), pch=".", main="Volcano plot",
     xlim=c(-4,4),xlab ="Log Fold-change",
     ylab = "-10*log10(FDR)",cex=1.5
     )
abline (h = -10*log10 (0.05), col = "blue")
abline (v = - log2 (2),col = "blue")
abline (v =  log2 (2),col = "blue")
text (-3, -10*log10 (0.05) -0.3 + 0.2,'FDR< 0.05',col="blue",cex=0.8)
grid()
#highlight only points
points(res.degs.up$lfc,-10*log10(res.degs.up$pval),labels="",pch="o",cex=0.5,col="red")
points(res.degs.down$lfc,-10*log10(res.degs.down$pval),labels="",pch="o",cex=0.5,col="green")
#highlight also gene names 
text(res.degs.up$lfc,-10*log10(res.degs.up$pval),labels=rownames(res.degs.up),cex=0.4,col="red",adj = c(0,0))
text(res.degs.down$lfc,-10*log10(res.degs.down$pval),labels=rownames(res.degs.down),cex=0.4 ,col="green",adj = c(1,1))

### Heatmap of DEGS 
meta.sub=meta[ meta$sample.type2 %in% c(group1, group3),]
column_ha = HeatmapAnnotation(Sample.type = meta.sub$sample.type2 )
Heatmap(exp1.degs, row_names_gp = gpar(fontsize = 2.5),name="Exp", column_names_gp = gpar(fontsize = 6), top_annotation = column_ha)
Heatmap(t(scale(t(exp1.degs))),row_names_gp = gpar(fontsize = 2.5),name="Z-score", column_names_gp = gpar(fontsize = 6), top_annotation = column_ha)




