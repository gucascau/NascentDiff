

### The script is to use poisson test to check the CRISPR and generate the viocano plot as well as the 

library(fitdistrplus) 
library(clusterProfiler)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")

library(org.Mm.eg.db)
keytypes(org.Mm.eg.db) 


library(ggplot2)
library(ggrepel)

### using their generated results and draw the volcano plot and GO enrichment

setwd("/Users/xinwang/Documents/Projects/RNA_DDR/Jass_Zulong/Results/Updated_0818_PossionTest/CRISPERscreen")

### Use the raw reads
RawReads<-read.table("062719_Nor_Raw_ZC.txt", header = T, row.names = 1)
attach(RawReads)

dim(RawReads)
## Generate the poisson test

# run the poisson test for two samples
mean(RawReads$high1_S3.nor/RawReads$unsort1_S1.nor)

Test<-poisson.test(x=4,T=5,alternative = "greater",conf.level = 0.95)
Test$p.value

head(RawReads)

pvalue=apply(RawReads,1,function(x) poisson.test(x =x[i,4], T=x[i,8], alternative = "greater",conf.level = 0.95)$p.value)
Pvalue<-NULL
for (i in 1:nrow(RawReads)){
  y=poisson.test(x = round(max(RawReads[i,4],RawReads[i,8])), r=min(RawReads[i,4],RawReads[i,8]) , alternative = "greater")
  
  Pvalue[i]<-y$p.value
}
?p.adjust
Padjust<-p.adjust(Pvalue,method = "fdr" )
Padjust
?p.adjust

LogFC<-log2(RawReads$high1_S3.nor/RawReads$unsort1_S1.nor)

table<-cbind(RawReads,Pvalue,Padjust,LogFC)

### draw the volcano plot

Proteinneddylation<-  c("Nae1","Rbx1","Rnf7","Uba3","Ube2m")

ggplot(table,aes(x=LogFC,y=-log10(Padjust)))+
  geom_point(size=0.5,color= ifelse(Padjust<=0.01 & LogFC>=0.58,"red",
                                    ifelse(Padjust<=0.01 & LogFC<=-0.58,"blue","gray")))+ 
  geom_text_repel(aes(label=ifelse(Padjust<=0.0000001 ,as.character(table$gene),'')),segment.size=0.01,
             hjust=0.2, size=3)+
  
  theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +ylim(0,30) + scale_x_continuous(limits=c(-6,6),breaks =c(-6,-4,-2,0,2,4,6))

### with a highlighed genes in protein neddylation
Proteinneddylation<-  c("Nae1","Rbx1","Rnf7","Uba3","Ube2m")

Proteinddeneddylation<- c(Cops3/Cops4/Cops7b/Cops8/Gps1)
table$genelabels <- ""
table$genelabels <- ifelse(table$gene == "Nae1"
                                  | table$gene == "Rbx1"
                                  | table$gene == "Rnf7"
                                  | table$gene == "Uba3"
                                  | table$gene == "Ube2m"
                                  | table$gene == "Ube2f"
                                  | table$gene == "Atm"
                                  | table$gene == "Nedd8"
                                  | table$gene =="Cul4b",
                                #  | table$gene == "Cops3"
                                #  | table$gene == "Cops4"
                                #  | table$gene == "Cops7b"
                                #  | table$gene == "Cops8" 
                                 # | table$gene == "Gps1", 
                           T, F)
table$genelabels 

ggplot(table,aes(x=LogFC,y=-log10(Padjust)))+
  geom_point(size=0.5,color= ifelse(Padjust<=0.01 & LogFC>=0.58,"red",
                                    ifelse(Padjust<=0.01 & LogFC<=-0.58,"blue","gray")))+ 
  #geom_text_repel(aes(label=ifelse(Padjust<=0.0000001 ,as.character(table$gene),'')),segment.size=0.01,
  #                hjust=0.2, size=3)+
  geom_text_repel(aes(label = ifelse(genelabels == TRUE, as.character(gene),"")),segment.size=0.01, hjust=0.2, size=3) + 
  
  
  theme_bw()+
  theme(plot.background = element_blank(),
        text = element_text(size=16),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + scale_x_continuous(limits=c(-6,6),breaks =c(-6,-4,-2,0,2,4,6)) +ylim(0,20)



BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

#https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
#res1 <- read.table("your_file_with_gene_fc_pvalue",header = TRUE,row.names = 1,sep = "\t")



AUDEG<- table[table$Padjust <= 0.01 &  table$LogFC >=0.58, ]

ADDEG<- table[table$Padjust <= 0.01 &  table$LogFC<=-0.58, ]

nrow(AUDEG)
nrow(ADDEG)

length(unique(AUDEG$gene))
length(unique(ADDEG$gene))

AgoUp<-enrichGO(unique(AUDEG$gene), OrgDb = "org.Mm.eg.db", ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.5,keyType = 'SYMBOL')

AgoUPSimply<-simplify(AgoUp)

AgoDown<-enrichGO(unique(ADDEG$gene), OrgDb = "org.Mm.eg.db", ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.5,keyType = 'SYMBOL')
AgoDownSimply<-simplify(AgoDown)

### Generate final table for the expression with P value P adjust and FC

write.csv(table, file="CRISPR_CompareHighToUnsorted.csv")

### Generate final Upregulated gene ontology and Down-regulated gene ontology

write.csv(AgoUp, file="ActivatedGenesGOenrichement_CompareHighToUnsorted.csv")

write.csv(AgoDown, file="RepressedGenesGOenrichement_CompareHighToUnsorted.csv")

### simplifed GO enrichment

write.csv(AgoUPSimply, file="ActivatedGenesSimplifiedGOenrichement_CompareHighToUnsorted.csv")

write.csv(AgoDownSimply, file="RepressedGenesSimplifiedGOenrichement_CompareHighToUnsorted.csv")

















simplify(AgoDown)
#AUDEG<- Genelist[Genelist$pos.p.value <= 0.05 &  Genelist$pos.lfc <=0.6667, ]


### check whether they follows the poisson test


library(dplyr)
library(purrr)

df %>% mutate(rate = map2(count, time, poisson.test))

lambda.est <- mean(RawReads$high1_S3.nor)

## check whether following poisson distribution
descdist(data=RawReads$high1_S3.nor,discrete=FALSE)

### fit the distirbution with pois
mean(RawReads$high1_S3.nor)
HighNumber<-RawReads$high1_S3.nor
fit_poi<-fitdist(HighNumber,"gamma",method="mme")
summary(fit_poi)
?fitdist
plotdist(HighNumber,histo = TRUE, demp=TRUE)

summary(fit_poi)
## 
gofstat(fit_poi)



#AUDEG<- Genelist[Genelist$pos.p.value <= 0.05 &  Genelist$pos.lfc >=1.5, ]

#AUDEG<- Genelist[Genelist$pos.p.value <= 0.05 &  Genelist$pos.lfc <=0.6667, ]
#AUDEG<- res$table[res$table$PValue <= 0.05 &  res$table$logFC >=0.58, ]
#ADDEG<- res$table[res$table$PValue <= 0.05 &  res$table$logFC <=- 0.58, ]
nrow(AUDEG)
nrow(ADDEG)

ATCO<-merge(out,y2, by="row.names")

write.table(ATCO,file="ZCCombinedSelf_All_v0228_v0228_cpm.txt")
?res
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

write.table(ATCO,file="ZCCombinedSelf_All_v0228_v0228_withoutLibrary_TMM.txt")

