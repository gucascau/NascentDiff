
### Generate the heatmap using the Z-score
library(ggplot2)
install.packages("pheatmap")


### here is for the latest histone genes:

setwd("/Users/xinwang/Documents/Projects/RNA_DDR/Jass_Zulong/Results/v0721_norDNA")
df <- read.table(file = "HistoneGeneExpressionNormalizedByLength_ForheatMapFig6F.txt", sep = "\t", header = T, row.names = 1,na.strings = NA)
attach(df)

library(pheatmap)
library(gplots)

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

