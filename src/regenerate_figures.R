# regenerate Figure 6 and 7

# Fig. 6B
# barplot
cpm <- read.csv("read_count_normalized_by_ERCC_total.csv", row.names = 1)

cpm.total <- colSums(cpm[,2:5])
cpm.rRNA <- colSums(cpm[grep("ribosomal", cpm$Gene), 2:5])
cpm.norRNA <- colSums(cpm[-grep("ribosomal", cpm$Gene), 2:5])
cpm.protein <- colSums(cpm[cpm$Gene %in% protein_coding_genes, 2:5])

e <- as.numeric(cpm.total)
e.mean <- c(mean(e[1:2]), mean(e[3:4])) *1e-6
e.sd <- c(sd(e[1:2]), sd(e[3:4])) * 1e-6

e1 <- as.numeric(cpm.rRNA)
e1.mean <- c(mean(e1[1:2]), mean(e1[3:4])) *1e-6
e1.sd <- c(sd(e1[1:2]), sd(e1[3:4])) * 1e-6

e2 <- as.numeric(cpm.protein)
e2.mean <- c(mean(e2[1:2]), mean(e2[3:4])) *1e-6
e2.sd <- c(sd(e2[1:2]), sd(e2[3:4])) * 1e-6

mat <- as.matrix(cbind(e.mean, e1.mean, e2.mean))
mat.sd <- cbind(e.sd, e1.sd, e2.sd)

dir.create("Figures/")
pdf("Figures/Figure6B.pdf")
par(mar=c(2,5,1,2))
bp <- barplot(mat, axes = F, beside = T, ylim = c(0,100), names.arg = rep("",3), 
              legend.text = FALSE, space = c(0.2, 0.8), col = c("blue", "red"))
points(rbind(bp - 0.2, bp + 0.2)[c(1,3,2,4), ], cbind(e, e1, e2) * 1e-6, pch = 19)
axis(1, at=c(0,colMeans(bp), colMeans(bp)[3] + 2), labels=F, lwd=3, pos=0)
mtext(text=c("","","","","",""), at=c(0,bp,bp[length(bp)] + 2), line=3, side=1, cex=3)
#axis(2, at = c(0,10,20,30,40,50), cex.axis=2, las=2, lwd=3, labels = T)
axis(2, cex.axis=2, las=2, lwd=3, labels = T)
arrows(bp, mat - mat.sd, bp, mat + mat.sd, angle=90, code=3, length = 0.05, lwd=3)
dev.off()

# changes in percentages
(mat[1,1] - mat[2,1]) / mat[1,1]
(mat[1,2] - mat[2,2]) / mat[1,2]
(mat[2,3] - mat[1,3]) / mat[1,3]

# poisson test
poisson.test(x = round(max(c(mean(e[1:2]), mean(e[3:4])))), 
             r = min(c(mean(e[1:2]), mean(e[3:4]))), 
             alternative = "greater")$p.value
poisson.test(x = round(max(c(mean(e1[1:2]), mean(e1[3:4])))), 
             r = min(c(mean(e1[1:2]), mean(e1[3:4]))), 
             alternative = "greater")$p.value
poisson.test(x = round(max(c(mean(e2[1:2]), mean(e2[3:4])))), 
             r = min(c(mean(e2[1:2]), mean(e2[3:4]))), 
             alternative = "greater")$p.value

# t-test (log2)
t.test(log2(e[1:2]), log2(e[3:4]), alternative = "greater")$p.value
t.test(log2(e1[1:2]), log2(e1[3:4]), alternative = "greater")$p.value
t.test(log2(e2[1:2]), log2(e2[3:4]), alternative = "less")$p.value

# Figure 6C

cpm.rRNA.sep <- cpm[grep("ribosomal", cpm$Gene), ]
cpm.rRNA.sep <- cpm.rRNA.sep[c(2,1,3), ]
rownames(cpm.rRNA.sep) <- cpm.rRNA.sep$Gene
cpm.rRNA.sep <- cpm.rRNA.sep[, 2:5]

e <- as.numeric(cpm.rRNA.sep[1,])
e.mean <- c(mean(e[1:2]), mean(e[3:4])) *1e-6
e.sd <- c(sd(e[1:2]), sd(e[3:4])) * 1e-6

e1 <- as.numeric(cpm.rRNA.sep[2,])
e1.mean <- c(mean(e1[1:2]), mean(e1[3:4])) *1e-6
e1.sd <- c(sd(e1[1:2]), sd(e1[3:4])) * 1e-6

e2 <- as.numeric(cpm.rRNA.sep[3,])
e2.mean <- c(mean(e2[1:2]), mean(e2[3:4])) *1e-6
e2.sd <- c(sd(e2[1:2]), sd(e2[3:4])) * 1e-6

emat <- as.matrix(cbind(e.mean, e1.mean, e2.mean))
mat.sd <- cbind(e.sd, e1.sd, e2.sd)


pdf("Figures/Figure6C.pdf", height = 4, width = 6)
par(mar=c(2,5,1,2))
bp <- barplot(mat, axes = F, beside = T, ylim = c(0,40), names.arg = rep("",3), 
              legend.text = FALSE, space = c(0.2, 0.8), col = c("blue", "red"))
points(rbind(bp - 0.2, bp + 0.2)[c(1,3,2,4), ], cbind(e, e1, e2) * 1e-6, pch = 19)
axis(1, at=c(0,colMeans(bp), colMeans(bp)[3] + 2), labels=F, lwd=3, pos=0)
mtext(text=c("","","","","",""), at=c(0,bp,bp[length(bp)] + 2), line=3, side=1, cex=3)
#axis(2, at = c(0,10,20,30,40,50), cex.axis=2, las=2, lwd=3, labels = T)
axis(2, cex.axis=2, las=2, lwd=3, labels = T)
arrows(bp, mat - mat.sd, bp, mat + mat.sd, angle=90, code=3, length = 0.05, lwd=3)
dev.off()

# changes in percentages
(mat[1,1] - mat[2,1]) / mat[1,1]
(mat[1,2] - mat[2,2]) / mat[1,2]
(mat[1,3] - mat[2,3]) / mat[1,3]

# poisson test
poisson.test(x = round(max(c(mean(e[1:2]), mean(e[3:4])))), 
             r = min(c(mean(e[1:2]), mean(e[3:4]))), 
             alternative = "greater")$p.value
poisson.test(x = round(max(c(mean(e1[1:2]), mean(e1[3:4])))), 
             r = min(c(mean(e1[1:2]), mean(e1[3:4]))), 
             alternative = "greater")$p.value
poisson.test(x = round(max(c(mean(e2[1:2]), mean(e2[3:4])))), 
             r = min(c(mean(e2[1:2]), mean(e2[3:4]))), 
             alternative = "greater")$p.value

# t test (log2)
t.test(log2(e[1:2]), log2(e[3:4]), alternative = "greater")$p.value
t.test(log2(e1[1:2]), log2(e1[3:4]), alternative = "greater")$p.value
t.test(log2(e2[1:2]), log2(e2[3:4]), alternative = "greater")$p.value


# Figure 6E

# protein coding
cpm.mean.pro <- cpm.mean[rownames(cpm.mean) %in% protein_coding_genes, ]

pdf("Figures/Figure6D.pdf")
plot(cpm.mean.pro[, "log2FC"], -log10(as.numeric(cpm.mean.pro[, "adj.P"])), col = cpm.mean.pro[, "color"], 
     pch = 19, cex = 0.3, ylim = c(0,20), xlim = c(-8, 8))
dev.off()
length(which(cpm.mean.pro[, "color"] == "red"))
length(which(cpm.mean.pro[, "color"] == "blue"))


# Figure 6F
library(pheatmap)
cols = colorRampPalette(c("blue","white","red"))(50)
DEG_table <- rbind(exp_table[exp_table$ENS %in% up_protein$ENS, c("noIR1","IR1","noIR2","IR2")],
                   exp_table[exp_table$ENS %in% down_protein$ENS, c("noIR1","IR1","noIR2","IR2")])
#rm <- which(rowSums(DEG_table) == 0)
pdf("Figures/Figure6F.pdf")
pheatmap(as.matrix(log2(DEG_table + 0.001)), scale = "row",  color = cols,
         fontsize_col=15, cluster_rows=FALSE, cluster_cols = FALSE,treeheight_row = 0, treeheight_col = 0,
         show_rownames = F, show_colnames = T)
dev.off()

# Figure 6G

# read in detailed table from Xin
exp_table <- read.table("read_count_normalized_by_ERCC_total_poisson_test_detailed.txt", header = T, sep = "\t")

up_protein <- exp_table[exp_table$Type == "protein_coding" & exp_table$log2FC > 0.58 & exp_table$adj.P < 0.05, ]
down_protein <- exp_table[exp_table$Type == "protein_coding" & exp_table$log2FC < -0.58 & exp_table$adj.P < 0.05, ]

write.csv(up_protein, file = "up-regulated_protein_coding_gene.csv")
write.csv(down_protein, file = "down-regulated_protein_coding_gene.csv")

library(clusterProfiler)
library(org.Mm.eg.db)
up_protein_BP <- enrichGO(gene = up_protein$ENS, OrgDb = "org.Mm.eg.db", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                          ont = "BP", qvalueCutoff = 0.5,keyType = "ENSEMBL")
#up_protein_BP_sim <- simplify(up_protein_BP)
up_protein_BP_sim <- setReadable(up_protein_BP_sim, OrgDb = "org.Mm.eg.db")
write.csv(data.frame(up_protein_BP_sim), file = "up_regulated_protein_GOBP_enrichment.csv")

down_protein_BP <- enrichGO(gene = down_protein$ENS, OrgDb = "org.Mm.eg.db", pvalueCutoff = 0.05, pAdjustMethod = "BH",
                          ont = "BP", qvalueCutoff = 0.5,keyType = "ENSEMBL")
#down_protein_BP_sim <- simplify(down_protein_BP)
down_protein_BP_sim <- setReadable(down_protein_BP_sim, OrgDb = "org.Mm.eg.db")
write.csv(data.frame(down_protein_BP_sim), file = "down_regulated_protein_GOBP_enrichment.csv")

# Figure 7A
library(xlsx)
normReadCount <- read.xlsx("NormalizedReadCounts_Figure7A_Xin.xlsx", sheetIndex = 1)

poisson.test(x = round(max(c(mean(as.numeric(normReadCount[1,2:3]), mean(as.numeric(normReadCount[1,4:5])))))), 
             r = min(c(mean(as.numeric(normReadCount[1,2:3])), mean(as.numeric(normReadCount[1,4:5])))), 
             alternative = "greater")$p.value
poisson.test(x = round(max(c(mean(as.numeric(normReadCount[2,2:3]))), mean(as.numeric(normReadCount[2,4:5])))), 
             r = min(c(mean(as.numeric(normReadCount[2,2:3])), mean(as.numeric(normReadCount[2,4:5])))), 
             alternative = "greater")$p.value
poisson.test(x = round(max(c(mean(as.numeric(normReadCount[3,2:3]))), mean(as.numeric(normReadCount[3,4:5])))), 
             r = min(c(mean(as.numeric(normReadCount[3,2:3])), mean(as.numeric(normReadCount[3,4:5])))), 
             alternative = "greater")$p.value
