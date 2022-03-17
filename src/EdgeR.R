#Author:Xin Wang 
#email: xin.wang@childrens.harvard.edu
#PI: Kaifu Chen

### Meanwhile, when we used the EdgeR, no more further normalization was used to detect the differential expression genes.
library(edgeR)
library(clusterProfiler)
library(org.Mm.eg.db)
keytypes(org.Mm.eg.db) 

setwd("~/Documents/Projects/RNA_DDR/Jass_Zulong/Results/v0721_norDNA/")

### For all the datasets, we did include the rDNA
mydata_norDNA<-read.table("ZCCombined2_noirDNA_selfNormForEdgeR2.txt", header=T,sep = "\t",row.names=1)

### quantile normalize the dataset. this step was igonored 
#counts<-normalizeQuantiles(mydata_norDNA)

counts<-mydata_norDNA

### set up the factor of each sample
group_all<- as.factor(c("Ctl", "Trt","Ctl","Trt"))
#x <- as.factor(c("Ctl", "Trt","Ctl","Trt"))
group_all
y<-DGEList(counts=counts, group=group_all)

### filter low counts genes
#keep<-rowSums(y$counts) >50

### Genes that have very low counts across all the libraries should be removed prior to downstream analysis. This is justified on both biological and statistical grounds. From biological point of view, a gene must be expressed at some minimal level before it is likely to be translated into a protein or to be considered biologically important. From a statistical point of view, genes with consistently low counts are very unlikely be assessed as significantly DE because low counts do not provide enough statistical evidence for a reliable judgement to be made. Such genes can therefore be removed from the analysis without any loss of information.

#####As a rule of thumb, we require that a gene have a count of at least 10–15 in at least some libraries before it is considered to be expressed in the study. We could explicitly select for genes that have at least a couple of counts of 10 or more, but it is slightly better to base the filtering on count-per-million (CPM) values so as to avoid favoring genes that are expressed in larger libraries over those expressed in smaller libraries. For the current analysis, we keep genes that have CPM values above 0.5 in at least two libraries:
# use cpm normalization
keep<-rowSums(cpm(y)>0.5) >2

table(keep)

### DGElList object is subsetted toreatin only the no-filtered genes.
genelist.filted<-y[keep,,keep.lib.sizes=FALSE]

##### Normalization for composition bias, we used the "none " method that considered normalization factor as 1.  

genelist.norm<-calcNormFactors(genelist.filted,method="none")

#### In this way, we did not considered the library size and norm factors for the dataset. 
genelist.norm$samples

#genelist.norm.log<-cpm(genelist.norm, log=TRUE)


### The follow is to Explore differences between libraries
label_name<-c("Ctl", "Trt","Ctl","Trt")

pch <- c(rep(21,4),rep(22,4))
plotMDS(genelist.norm, labels = label_name, pch = 22, bg = par("bg"),cex = 1, col=c(rep("black",4), rep("red",4)))


?plotMDS


mds <- plotMDS(genelist.norm)
toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, Group = c("Ctl", "Trt","Ctl","Trt"))
library(ggplot2)
ggplot(toplot, aes(Dim1, Dim2,col=c(rep("black",4), rep("red",4))))+ geom_point() +
  #geom_text( label = label_name, nudge_x = 0.1, nudge_y = 0.1,  check_overlap = T)+
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

##### check the batch effects
#group_all
#$logCPM_Batch <- removeBatchEffect(genelist.norm.log, batch=group_all)


#plotMDS(logCPM_Batch)
#?plotMDS

####### Pre-treatment



#group_list<-factor(c(rep("Ctret",3), rep("Atret",3), rep("Btret",3), rep("Con",3)))
#y <- DGEList(counts=counts, group=group_list)


#### design matrix ：
design <- model.matrix(~0+group_all)

colnames(design) <-levels(group_all)

design

##### 估计离散值 （Dispersion） Dispersion estimation
genelist.Disp <- estimateDisp(genelist.norm, design, robust = TRUE)

plotBCV(genelist.Disp)

fit <- glmQLFit(genelist.Disp, design, robust=TRUE)

head(fit$coefficients)


### Testing for differential expression
cntrvsKD <- makeContrasts(Trt-Ctl, levels=design)
res <- glmQLFTest(fit, contrast=cntrvsKD)

Tablle_new<-res$fitted.values

plotMD(res,ylim=c(-6,6)) 
?plotMD
ig.edger <- res$table[p.adjust(res$table$PValue, method = "BH") < 1, ]

#### Generate the pvalue, pajust value, and logFC, The top DE genes can be viewed with topTags:

out<- topTags(res,  adjust.method = "BH", sort.by = "PValue", n="Inf")$table


### identify the upregualted genes and downregualted genes

AUDEG<- res$table[p.adjust(res$table$PValue, method = "BH") <= 0.05 &  res$table$logFC >=0.58, ]
ADDEG<- res$table[p.adjust(res$table$PValue, method = "BH") <= 0.05 &  res$table$logFC <=-0.58, ]
#AUDEG<- res$table[res$table$PValue <= 0.05 &  res$table$logFC >=0.58, ]
#ADDEG<- res$table[res$table$PValue <= 0.05 &  res$table$logFC <=- 0.58, ]
nrow(AUDEG)
nrow(ADDEG)


### Then dara the heatmap clustering


write.table(AUDEG, file = "ZCCombinedSelf_Upregulated_v0228.txt")

write.table(AUDEG, file = "ZCCombinedSelf_Downregulated_v0228.txt")

AgoDEG<-enrichGO(row.names(AUDEG), OrgDb = "org.Mm.eg.db", ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.5,keyType = 'ENSEMBL')

write.table(AgoDEG,file = "ZCCombinedSelf_UpregulatedGO_v0228.txt")
write.table(ADDEG,file="ZCCombinedSelf_Downregulated_v0228.txt")

#write.table(AUDEG,file="Jun_OverExpress_UP_Gens_two.txt")
### GO enrichment of Upupregulted_genes and down-regulated genes 
AgoUp<-enrichGO(row.names(AUDEG), OrgDb = "org.Mm.eg.db", ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.5,keyType = 'ENSEMBL')
simplify(AgoUp)
AgoDown<-enrichGO(row.names(ADDEG), OrgDb = "org.Mm.eg.db", ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.5,keyType = 'ENSEMBL')
head(AgoUp)
head(AgoDown)

write.table(AgoUp, file="ZCCombinedSelf_UpregulatedGO_v0228.txt")
write.table(AgoDown, file="ZCCombinedSelf_DownregulatedGO_v0228.txt")

ATCO<-merge(out,res$fitted.values, by="row.names")

write.table(ATCO,file="ZCCombinedSelf_All_v0228_v0228.txt")




### Here is for heatmap:
dev.off()
## The original heatmap --- one
logCPM <- cpm(counts, prior.count=2, log=TRUE)

rownames(logCPM) <- counts$genes$Symbol 

colnames(logCPM)
#colnames(logCPM) <- paste(counts$samples$group, 1:2, sep="-")

## Use the EdgR heatmap --- two

logCPM <- cpm(genelist.Disp, prior.count=2, log=TRUE)

rownames(logCPM) <- genelist.Disp$genes$Symbol 
colnames(logCPM) <- paste(genelist.Disp$samples$group, 1:2, sep="-")

### we shoould not 
logCPM <- t(scale(t(logCPM)))

nrow(logCPM)
library(gplots)
col.pan <- colorpanel(100, "blue", "white", "red")

heatmap.2(as.matrix(genelist.Disp), scale = "row", col = bluered(1000), Rowv = TRUE, Colv = FALSE,
          #Rowv = Rowv, Colv = Colv,Rowv = TRUE
          trace = "none", density.info = "none")


heatmap.2(logCPM, col=col.pan,  scale="row",
             trace="none", dendrogram="none", cexRow=1, cexCol=1.4, density.info="none", Rowv = TRUE, Colv = FALSE,
            margin=c(10,9), lhei=c(2,10), lwid=c(2,6))
