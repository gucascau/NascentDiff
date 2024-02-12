

##### draw the figure with the histone gene expression


### Generate the heatmap using the Z-score
library(ggplot2)

library(pheatmap)
library(gplots)
library(ggrepel)
#install.packages("pheatmap")


### here is for the latest histone genes:

setwd("/Users/xinwang/Documents/Projects/RNA_DDR/Jass_Zulong/Results/Updated_0818_PossionTest/ERCC_normalized_read_count_and_differentially_expressed_genes/ProcessingResults")
df <- read.table(file = "read_count_normalized_by_ERCC_total_poisson_test_detailed.txt", sep = "\t", header = T, row.names = 1,na.strings = NA)
attach(df)

### extract all the protein coding genes

Proteincoding<-df[df$Type=="protein_coding",]

## then draw the figure 6D scatterplot ranked by the mean non irradiation value

### we generated the RPKM value


RPKM<-log2((Proteincoding$noIR1.nl+Proteincoding$noIR2.nl)/2)

FinalCoding<-cbind(Proteincoding,RPKM)
attach(FinalCoding)

#pdf(file = "ScatterPlot_MeanExpression_WT.pdf")
ggplot(FinalCoding,aes(x=RPKM,y=log2FC))+
  geom_point(size=1,color= ifelse(log2FC>0.58,"red",
                                  ifelse(log2FC<= -0.58,"blue","gray")))+ 
  
  
 # geom_text_repel(aes(label=ifelse((RPKM>0.5 & (log2FC >=1 | log2FC<=-1)),as.character(GeneID),'')),segment.size=0.01,
 #                 hjust=0.2, size=3,max.overlaps =25)+
  geom_hline(yintercept=0.58, linetype="dashed", color = "black")+
  geom_hline(yintercept=-0.58, linetype="dashed", color = "black")+
  
  geom_hline(yintercept=0,  color = "black")+
  theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())  +xlab("Mean Gene Expression") + ylab("log2(Fold Change) (IR vs Control)")+
scale_x_continuous(breaks=seq(-16, 8, 4),limits = c(-16,8)) +scale_y_continuous(breaks=seq(-6, 6, 3),limits = c(-6,6))

dev.off()


## We generated figure 6F  Histone gene expression


### only select the row with histone genes
HistoneGenes<-df[ grep("Hist",df$GeneID,),]

### setup the expression pattern of these histone genes
HistoneGeneExp<-df[ grep("Hist",df$GeneID,),c(15,17,16,18)]

### setup the new name of these histone genes
rownames(HistoneGeneExp)<-df[ grep("Hist",df$GeneID,),c(2)]

attach(HistoneGeneExp)

### 
### sorted the histone genes by expression
new_mtx <- HistoneGeneExp [which(rowSums(HistoneGeneExp) > 0), ]  # works much faster than below
new_mtx_sorted <-new_mtx[order(new_mtx[,1],decreasing=TRUE),]
#new_mtx  <- apply (df, 1, function(x) any(x==0) )

#new_mtx <-filter_if(df, is.numeric, all_vars((.) != 0)) 

cols=bluered(100)
#cols = colorRampPalette(c("blue", "red"))(100)

#pheatmap(as.matrix(new_mtx), color=cols, show_rownames = F,cluster_cols=T,cluster_rows=F, scale="row", cex=1,clustering_distance_rows = "euclidean", cex=1, clustering_distance_cols = "euclidean", clustering_method = "complete",border_color = FALSE)
pheatmap(as.matrix(log2(new_mtx_sorted+1)), scale = "none", clustering_distance_row = "correlation", color = cols,
              clustering_method = "complete",
              fontsize_col=15, cluster_rows=F, cluster_cols = F,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = T)


dev.off()
### We generated figure 6B the top 100 genes


### identify differential expressed genes

FinalCodingSorted<-FinalCoding[order(FinalCoding$RPKM,FinalCoding$IR1.nl,decreasing=TRUE),]
DEGene<-FinalCodingSorted[FinalCodingSorted$adj.P<=0.05 & (FinalCodingSorted$log2FC >=0.58 | FinalCodingSorted$log2FC <=-0.58),]

### sorted the top 50 genes sorted by FinalCoding RPKM


DEGeneSorted <-DEGene[order(DEGene$RPKM,DEGene$IR1.nl,decreasing=TRUE),]
### print the DEGGensSorted into a file
write.table(DEGeneSorted,file="DEGeneListSorted.txt",sep = "\t")

TopDEGeneSortedAinf<-DEGeneSorted[1:100,]
TopDEGeneSorted<-DEGeneSorted[1:100,c(15,17,16,18)]
rownames(TopDEGeneSorted)<-DEGeneSorted[1:100,c(2)]

### we used the expressed z-score to show heatmap
### figure B
pheatmap(as.matrix(log2(TopDEGeneSorted)), scale = "row",  color = cols,
              fontsize_col=15, cluster_rows=FALSE, cluster_cols = FALSE,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = T)

pheatmap(as.matrix(log2(TopDEGeneSorted)), scale = "none", clustering_distance_row = "correlation", color = cols,
         clustering_method = "complete",
         fontsize_col=15, cluster_rows=F, cluster_cols = F,treeheight_row = 0, treeheight_col = 0,
         show_rownames = T, show_colnames = T)

### generate the final size distribution (Barplot) Fig 6C

## UpRegulated gene size
UpRegulated<- FinalCodingSorted[(FinalCodingSorted$adj.P<=0.05 & FinalCodingSorted$log2FC >=0.58 ),]
NumUpRegulated<-nrow(UpRegulated)
NumUpRegulated
### The gene size of UpRegulated 
UpRegulatedSize<-cbind(replicate(NumUpRegulated,"Activated"),UpRegulated$GeneLength)

## DownRegulated gene size
DownRegulated<- FinalCodingSorted[(FinalCodingSorted$adj.P<=0.05 & FinalCodingSorted$log2FC <= -0.58 ),]
NumDownRegulated<-nrow(DownRegulated)
NumDownRegulated
DownRegulatedSize<-cbind(replicate(NumDownRegulated,"Repressed"),DownRegulated$GeneLength)

## Gene without any change
Nochange<-FinalCodingSorted[-FinalCodingSorted$adj.P<=0.05 & (FinalCodingSorted$log2FC >=0.58 | FinalCodingSorted$log2FC <=-0.58),]
NumNochange<-nrow(Nochange)
NochangeSize<-cbind(replicate(NumNochange,"nochange"),Nochange$GeneLength)


## Top 100 expressed genes

## Up 
UPTop100RegulatedGene<-TopDEGeneSortedAinf[(TopDEGeneSortedAinf$adj.P<=0.05 & TopDEGeneSortedAinf$log2FC >=0.58 ),]
# the gene size of Top100 Up regulated genes
UPTop100RegulatedGeneSize<-cbind(replicate(nrow(UPTop100RegulatedGene),"Top100Activated"),UPTop100RegulatedGene$GeneLength)

## Down 
DownTop100RegulatedGene<-TopDEGeneSortedAinf[(TopDEGeneSortedAinf$adj.P<=0.05 & TopDEGeneSortedAinf$log2FC <=-0.58 ),]
# the gene size of Top100 Up regulated genes
DownTop100RegulatedGeneSize<-cbind(replicate(nrow(DownTop100RegulatedGene),"Top100Repressed"),DownTop100RegulatedGene$GeneLength)

### then combine all the top 100 DEG, total DEGs
 
GenelengthDist<-rbind(UPTop100RegulatedGeneSize,DownTop100RegulatedGeneSize,UpRegulatedSize,DownRegulatedSize,NochangeSize)

write.table(GenelengthDist,file="GenelengthDist_differentCat.txt",sep = "\t")
colnames(GenelengthDist)<-c("Type","Length")
GenelengthDistDF<-as.data.frame(GenelengthDist)

### masure the length comparison Pvalue
DownTop100RegulatedGene$GeneLength
UPTop100RegulatedGene$GeneLength
? wilcox.test

#wilcox.test(sgs1$Frequency~sgs1$Type, data=sgs1,alternative = "less",var.equal = TRUE,paired=FALSE,correction=FALSE)
## comparing betwen top 100 repressed and acitvated DEG
wilcox.test(DownTop100RegulatedGene$GeneLength,UPTop100RegulatedGene$GeneLength,alternative = "less")

## comparing between top 100 repressed and total repressed
wilcox.test(DownTop100RegulatedGene$GeneLength,DownRegulated$GeneLength,alternative = "less")

### comparing between top 100 activated and total activated
wilcox.test(UPTop100RegulatedGene$GeneLength,UpRegulated$GeneLength,alternative = "less")


### comparing between repressed and non changed 
wilcox.test(DownRegulated$GeneLength,Nochange$GeneLength,alternative = "less")

### comparing between activated and no changed
wilcox.test(UpRegulated$GeneLength,Nochange$GeneLength,alternative = "greater")
## measure the total normalized read counts ranked by 

### generate the read distribution for DEG

# summary the top 100 for column four samples


### for the top 100 differential genes total normalized reads
Top100Reads<-colSums(DEGeneSorted[1:100,c(10,11,12,13)],na.rm = TRUE)
Top100Reads
### for the total DEG total normalized reads
DEGReads<-colSums(DEGeneSorted[,c(10,11,12,13)],na.rm = TRUE)
DEGReads
### for all normalized reads (protein coding)
TotalReads<-colSums(FinalCodingSorted[,c(10,11,12,13)],na.rm = TRUE)
TotalReads
### combine the three lists

CombineReadlist<-rbind(Top100Reads,TopDEGReads,TotalReads)
CombineReadlist



## for the top 100 normalized RPKM
colSums(DEGeneSorted[,c(15,17,16,18)],na.rm = TRUE)

























### This figure is for 

res$table[p.adjust(res$table$PValue, method = "BH") <= 0.05 &  res$table$logFC >=0.58, ]

### exclude the zero value
new_mtx <- df [which(rowSums(df) > 10), ]  # works much faster than below

#new_mtx  <- apply (df, 1, function(x) any(x==0) )

#new_mtx <-filter_if(df, is.numeric, all_vars((.) != 0)) 

cols=bluered(100)
#cols = colorRampPalette(c("blue", "red"))(100)

#pheatmap(as.matrix(new_mtx), color=cols, show_rownames = F,cluster_cols=T,cluster_rows=F, scale="row", cex=1,clustering_distance_rows = "euclidean", cex=1, clustering_distance_cols = "euclidean", clustering_method = "complete",border_color = FALSE)
out<-pheatmap(as.matrix(log2(new_mtx)), scale = "none", color = cols,
              clustering_method = "complete",
              fontsize_col=15, cluster_rows=F, cluster_cols = F,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = F)

geneclusters <- out[["kmeans"]][["cluster"]]

geneclusters[geneclusters==5]

library(gplots)
heatmap.2(as.matrix(df), scale = "row", col = bluered(1000), Rowv = TRUE, Colv = FALSE,
          #Rowv = Rowv, Colv = Colv,Rowv = TRUE
          trace = "none", density.info = "none")

?heatmap.2



### the following is to generate the heatmap for Cdk genes
df <- read.table(file = "Cyclin_dependentKinaseGeneForHeatMap_v0314.txt", sep = "\t", header = T, row.names = 1,na.strings = NA)
attach(df)

library(pheatmap)

### exclude the zero value
new_mtx <- df [which(rowSums(df) > 0), ]  # works much faster than below

#new_mtx  <- apply (df, 1, function(x) any(x==0) )

#new_mtx <-filter_if(df, is.numeric, all_vars((.) != 0)) 

cols=bluered(100)
#cols = colorRampPalette(c("blue", "red"))(100)

#pheatmap(as.matrix(new_mtx), color=cols, show_rownames = F,cluster_cols=T,cluster_rows=F, scale="row", cex=1,clustering_distance_rows = "euclidean", cex=1, clustering_distance_cols = "euclidean", clustering_method = "complete",border_color = FALSE)
out<-pheatmap(as.matrix(log2(new_mtx)), scale = "none", clustering_distance_row = "correlation", color = cols,
              clustering_method = "complete",
              fontsize_col=15, cluster_rows=F, cluster_cols = F,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = T)


### The followng is to generate a heatmap for  Negative regulation of cell cycle G1/S phase (83) 

df <- read.table(file = "GSEA/NegativeRegulationofCellCycleG1_Sphase_mostImportantGenesGeneExpForHeatMap.txt", sep = "\t", header = T, row.names = 1,na.strings = NA)
attach(df)

library(pheatmap)

### exclude the zero value
new_mtx <- df [which(rowSums(df) > 0), ]  # works much faster than below

#new_mtx  <- apply (df, 1, function(x) any(x==0) )

#new_mtx <-filter_if(df, is.numeric, all_vars((.) != 0)) 

cols=bluered(100)
#cols = colorRampPalette(c("blue", "red"))(100)

#pheatmap(as.matrix(new_mtx), color=cols, show_rownames = F,cluster_cols=T,cluster_rows=F, scale="row", cex=1,clustering_distance_rows = "euclidean", cex=1, clustering_distance_cols = "euclidean", clustering_method = "complete",border_color = FALSE)
out<-pheatmap(as.matrix(log2(new_mtx)), scale = "none", clustering_distance_row = "correlation", color = cols,
              clustering_method = "complete",
              fontsize_col=15, cluster_rows=F, cluster_cols = F,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = F)



### the following is to generate the heatmap for top 104 expressed
df <- read.table(file = "Top104Forheamap_v0228.txt", sep = "\t", header = T, row.names = 1,na.strings = NA)
attach(df)

library(pheatmap)

### exclude the zero value
new_mtx <- df [which(rowSums(df) > 10), ]  # works much faster than below

#new_mtx  <- apply (df, 1, function(x) any(x==0) )

#new_mtx <-filter_if(df, is.numeric, all_vars((.) != 0)) 

cols=bluered(100)
#cols = colorRampPalette(c("blue", "red"))(100)

#pheatmap(as.matrix(new_mtx), color=cols, show_rownames = F,cluster_cols=T,cluster_rows=F, scale="row", cex=1,clustering_distance_rows = "euclidean", cex=1, clustering_distance_cols = "euclidean", clustering_method = "complete",border_color = FALSE)
out<-pheatmap(as.matrix(log2(new_mtx)), scale = "row",  color = cols,
              fontsize_col=15, cluster_rows=FALSE, cluster_cols = FALSE,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = T)


### the following is to generate the heatmap for Cdk genes
df <- read.table(file = "Pol_v0228.txt", sep = "\t", header = T, row.names = 1,na.strings = NA)
attach(df)

library(pheatmap)

### exclude the zero value
new_mtx <- df [which(rowSums(df) > 10), ]  # works much faster than below

#new_mtx  <- apply (df, 1, function(x) any(x==0) )

#new_mtx <-filter_if(df, is.numeric, all_vars((.) != 0)) 

cols=bluered(100)
#cols = colorRampPalette(c("blue", "red"))(100)

#pheatmap(as.matrix(new_mtx), color=cols, show_rownames = F,cluster_cols=T,cluster_rows=F, scale="row", cex=1,clustering_distance_rows = "euclidean", cex=1, clustering_distance_cols = "euclidean", clustering_method = "complete",border_color = FALSE)
out<-pheatmap(as.matrix(log2(new_mtx)), scale = "none", clustering_distance_row = "correlation", color = cols,
              clustering_method = "complete",
              fontsize_col=15, cluster_rows=F, cluster_cols = F,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = T)

### exclude the zero value
new_mtx <- df [which(rowSums(df) > 10), ]  # works much faster than below

#new_mtx  <- apply (df, 1, function(x) any(x==0) )

#new_mtx <-filter_if(df, is.numeric, all_vars((.) != 0)) 

cols=bluered(100)
#cols = colorRampPalette(c("blue", "red"))(100)

#pheatmap(as.matrix(new_mtx), color=cols, show_rownames = F,cluster_cols=T,cluster_rows=F, scale="row", cex=1,clustering_distance_rows = "euclidean", cex=1, clustering_distance_cols = "euclidean", clustering_method = "complete",border_color = FALSE)
out<-pheatmap(as.matrix(log2(new_mtx)), scale = "none", clustering_distance_row = "correlation", color = cols,
              clustering_method = "complete",
              fontsize_col=15, cluster_rows=F, cluster_cols = F,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = T)


### The following is the heatmap for long noncoding DEG


### the following is to generate the heatmap for Cdk genes
df <- read.table(file = "TopLongNoncoding_v0228.txt", sep = "\t", header = T, row.names = 1,na.strings = NA)
attach(df)

library(pheatmap)

### exclude the zero value
new_mtx <- df [which(rowSums(df) > 10), ]  # works much faster than below

#new_mtx  <- apply (df, 1, function(x) any(x==0) )

#new_mtx <-filter_if(df, is.numeric, all_vars((.) != 0)) 

cols=bluered(100)
#cols = colorRampPalette(c("blue", "red"))(100)

#pheatmap(as.matrix(new_mtx), color=cols, show_rownames = F,cluster_cols=T,cluster_rows=F, scale="row", cex=1,clustering_distance_rows = "euclidean", cex=1, clustering_distance_cols = "euclidean", clustering_method = "complete",border_color = FALSE)
out<-pheatmap(as.matrix(log2(new_mtx)), scale = "row",  color = cols,
              fontsize_col=15, cluster_rows=FALSE, cluster_cols = FALSE,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = T)

out

### exclude the zero value
new_mtx <- df [which(rowSums(df) > 10), ]  # works much faster than below

#new_mtx  <- apply (df, 1, function(x) any(x==0) )

#new_mtx <-filter_if(df, is.numeric, all_vars((.) != 0)) 

cols=bluered(100)
#cols = colorRampPalette(c("blue", "red"))(100)

#pheatmap(as.matrix(new_mtx), color=cols, show_rownames = F,cluster_cols=T,cluster_rows=F, scale="row", cex=1,clustering_distance_rows = "euclidean", cex=1, clustering_distance_cols = "euclidean", clustering_method = "complete",border_color = FALSE)
out<-pheatmap(as.matrix(log2(new_mtx)), scale = "none", clustering_distance_row = "correlation", color = cols,
              clustering_method = "complete",
              fontsize_col=15, cluster_rows=F, cluster_cols = F,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = T)


### heatmap for Ribosomal proteins

df <- read.table(file = "RplGenes_v0306.txt", sep = "\t", header = T, row.names = 1,na.strings = NA)
attach(df)

library(pheatmap)

### exclude the zero value
new_mtx <- df [which(rowSums(df) > 10), ]  # works much faster than below

#new_mtx  <- apply (df, 1, function(x) any(x==0) )

#new_mtx <-filter_if(df, is.numeric, all_vars((.) != 0)) 

cols=bluered(100)
#cols = colorRampPalette(c("blue", "red"))(100)

#pheatmap(as.matrix(new_mtx), color=cols, show_rownames = F,cluster_cols=T,cluster_rows=F, scale="row", cex=1,clustering_distance_rows = "euclidean", cex=1, clustering_distance_cols = "euclidean", clustering_method = "complete",border_color = FALSE)
out<-pheatmap(as.matrix(log2(new_mtx)), scale = "row",  color = cols,
              fontsize_col=15, cluster_rows=T, cluster_cols = FALSE,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = T)

out

### exclude the zero value
new_mtx <- df [which(rowSums(df) > 10), ]  # works much faster than below

#new_mtx  <- apply (df, 1, function(x) any(x==0) )

#new_mtx <-filter_if(df, is.numeric, all_vars((.) != 0)) 

cols=bluered(100)
#cols = colorRampPalette(c("blue", "red"))(100)

#pheatmap(as.matrix(new_mtx), color=cols, show_rownames = F,cluster_cols=T,cluster_rows=F, scale="row", cex=1,clustering_distance_rows = "euclidean", cex=1, clustering_distance_cols = "euclidean", clustering_method = "complete",border_color = FALSE)
out<-pheatmap(as.matrix(log2(new_mtx)), scale = "none", clustering_distance_row = "correlation", color = cols,
              clustering_method = "complete",
              fontsize_col=15, cluster_rows=F, cluster_cols = F,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = T)


### heatmap for top5 normalized by gene length

df <- read.table(file = "Top5_normalizedByGenelength.txt", sep = "\t", header = T, row.names = 1,na.strings = NA)
attach(df)

library(pheatmap)

### exclude the zero value
new_mtx <- df [which(rowSums(df) > 10), ]  # works much faster than below

#new_mtx  <- apply (df, 1, function(x) any(x==0) )

#new_mtx <-filter_if(df, is.numeric, all_vars((.) != 0)) 

cols=bluered(100)
#cols = colorRampPalette(c("blue", "red"))(100)

#pheatmap(as.matrix(new_mtx), color=cols, show_rownames = F,cluster_cols=T,cluster_rows=F, scale="row", cex=1,clustering_distance_rows = "euclidean", cex=1, clustering_distance_cols = "euclidean", clustering_method = "complete",border_color = FALSE)
out<-pheatmap(as.matrix(log2(new_mtx)), scale = "row",  color = cols,
              fontsize_col=15, cluster_rows=T, cluster_cols = FALSE,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = T)

out

### exclude the zero value
new_mtx <- df [which(rowSums(df) > 10), ]  # works much faster than below

#new_mtx  <- apply (df, 1, function(x) any(x==0) )

#new_mtx <-filter_if(df, is.numeric, all_vars((.) != 0)) 

cols=bluered(100)
#cols = colorRampPalette(c("blue", "red"))(100)

#pheatmap(as.matrix(new_mtx), color=cols, show_rownames = F,cluster_cols=T,cluster_rows=F, scale="row", cex=1,clustering_distance_rows = "euclidean", cex=1, clustering_distance_cols = "euclidean", clustering_method = "complete",border_color = FALSE)
out<-pheatmap(as.matrix(log2(new_mtx)), scale = "row", clustering_distance_row = "correlation", color = cols,
              clustering_method = "complete",
              fontsize_col=15, cluster_rows=F, cluster_cols = F,treeheight_row = 0, treeheight_col = 0,
              show_rownames = T, show_colnames = T)

