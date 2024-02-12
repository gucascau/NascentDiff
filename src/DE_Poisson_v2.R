# normalized to mean of raw read counts
head(total_readCount)
colSums(total_readCount[,2:5])
mean_total_count <- mean(colSums(total_readCount[,2:5]))

# Method 1: normalized by ERCC total count
cpm <- total_readCount
cpm$noIR1 <- cpm$noIR1 / ERCC_sum["noIR1"]
cpm$noIR2 <- cpm$noIR2 / ERCC_sum["noIR2"]
cpm$IR1 <- cpm$IR1 / ERCC_sum["IR1"]
cpm$IR2 <- cpm$IR2 / ERCC_sum["IR2"]

scaleFactor <- mean_total_count / mean(colSums(cpm[,2:5]))

cpm$noIR1 <- cpm$noIR1 * scaleFactor
cpm$noIR2 <- cpm$noIR2 * scaleFactor
cpm$IR1 <- cpm$IR1 * scaleFactor
cpm$IR2 <- cpm$IR2 * scaleFactor

cpm.total <- colSums(cpm[,2:5])
cpm.rRNA <- colSums(cpm[grep("ribosomal", cpm$Gene), 2:5])
cpm.norRNA <- colSums(cpm[-grep("ribosomal", cpm$Gene), 2:5])
cpm.protein <- colSums(cpm[cpm$Gene %in% protein_coding_genes, 2:5])

write.csv(cpm, file = "read_count_normalized_by_ERCC_total.csv")

# Method 2: normalized by total read counts
cpm <- total_readCount
cpm$noIR1 <- cpm$noIR1 / colSums(cpm["noIR1"])
cpm$noIR2 <- cpm$noIR2 / colSums(cpm["noIR2"])
cpm$IR1 <- cpm$IR1 / colSums(cpm["IR1"])
cpm$IR2 <- cpm$IR2 / colSums(cpm["IR2"])

scaleFactor <- mean_total_count / mean(colSums(cpm[,2:5]))

cpm$noIR1 <- cpm$noIR1 * scaleFactor
cpm$noIR2 <- cpm$noIR2 * scaleFactor
cpm$IR1 <- cpm$IR1 * scaleFactor
cpm$IR2 <- cpm$IR2 * scaleFactor

cpm.total <- colSums(cpm[,2:5])
cpm.rRNA <- colSums(cpm[grep("ribosomal", cpm$Gene), 2:5])
cpm.norRNA <- colSums(cpm[-grep("ribosomal", cpm$Gene), 2:5])
cpm.protein <- colSums(cpm[cpm$Gene %in% protein_coding_genes, 2:5])

write.csv(cpm, file = "read_count_normalized_by_total_reads.csv")

# Method 3: normalized by median gene expression (based on total counts?)
cpm <- total_readCount
cpm$noIR1 <- cpm$noIR1 / median(cpm$noIR1)
cpm$noIR2 <- cpm$noIR2 / median(cpm$noIR2)
cpm$IR1 <- cpm$IR1 / median(cpm$IR1)
cpm$IR2 <- cpm$IR2 / median(cpm$IR2)

scaleFactor <- mean_total_count / mean(colSums(cpm[,2:5]))

cpm$noIR1 <- cpm$noIR1 * scaleFactor
cpm$noIR2 <- cpm$noIR2 * scaleFactor
cpm$IR1 <- cpm$IR1 * scaleFactor
cpm$IR2 <- cpm$IR2 * scaleFactor


# test
apply(cpm[, c("noIR1", "noIR2", "IR1", "IR2")], 2, median)
# Done!


cpm.total <- colSums(cpm[,2:5])
cpm.rRNA <- colSums(cpm[grep("ribosomal", cpm$Gene), 2:5])
cpm.norRNA <- colSums(cpm[-grep("ribosomal", cpm$Gene), 2:5])
cpm.protein <- colSums(cpm[cpm$Gene %in% protein_coding_genes, 2:5])

write.csv(cpm, file = "read_count_normalized_by_median.csv")

###################################################################
# Using Poisson Test to call differential gene expression
# Taking average of normalized read count
cpm.mean <- cbind(noIR_mean = rowMeans(cpm[, c("noIR1", "noIR2")]),
                  IR_mean = rowMeans(cpm[, c("IR1", "IR2")]))
rownames(cpm.mean) <- cpm$Gene

# compare average normalized read count between IR+ and IR-
# Taking the smaller value as the mean value in Possion Test, calculate P-value (one tail)
# remove genes with zero expression
rm <- which(rowMeans(cpm.mean[, c("noIR_mean", "IR_mean")]) == 0)
cpm.mean <- cpm.mean[-rm, ]
cpm.mean <- cbind(cpm.mean, 
                  log2FC = log2( (cpm.mean[,"IR_mean"] + 0.001) / (cpm.mean[, "noIR_mean"] + 0.001) ))

cpm.mean <- cbind(cpm.mean,
                  P = 0)
cpm.mean[, "P"] <-  apply(cpm.mean, 1, FUN = function(x) {
  poisson.test(x = round(max(x[c("noIR_mean","IR_mean")])), 
               r = min(x[c("noIR_mean","IR_mean")]), 
               alternative = "greater")$p.value
})


# Adjust P-value for all the genes

cpm.mean <- cbind(cpm.mean,
                  adj.P = 0)

cpm.mean[, "adj.P"] <- p.adjust(cpm.mean[,"P"], method = "fdr")

write.csv(cpm.mean, file = "read_count_normalized_by_ERCC_total_poisson_test.csv")
write.csv(cpm.mean, file = "read_count_normalized_by_total_reads_poisson_test.csv")
write.csv(cpm.mean, file = "read_count_normalized_by_median_poisson_test.csv")

# regenerate MA plot, volcano plot and barplot

# MA plot
maPlot(cpm.mean[, "noIR_mean"], cpm.mean[,"IR_mean"], 
       normalize = F, 
       ## only divides by the sum of counts in each sample and has nothing to do with the normalization factors
       pch = '.', 
       cex = 0.4, ylim = c(-5, 5),
       main = '',
)
grid(col = "blue")

# volcano plot
# volcano plot (protein-coding genes?)
cpm.mean <- as.data.frame(cpm.mean)
cpm.mean <- cbind(cpm.mean, 
                  color = "grey")
cpm.mean[cpm.mean$log2FC > 0.58 & cpm.mean$adj.P < 0.05, "color"] <- "red"
cpm.mean[cpm.mean$log2FC < -0.58 & cpm.mean$adj.P < 0.05, "color"] <- "blue"

upDEG <- cpm.mean[cpm.mean$log2FC > 0.58 & cpm.mean$adj.P < 0.05, ]
downDEG <- cpm.mean[cpm.mean$log2FC < -0.58 & cpm.mean$adj.P < 0.05, ]

write.csv(upDEG, file = "upDEG_IR_ERCC_total_poission_test.csv")
write.csv(downDEG, file = "downDEG_IR_ERCC_total_poission_test.csv")

# protein coding
cpm.mean.pro <- cpm.mean[rownames(cpm.mean) %in% protein_coding_genes, ]

pdf("Figures/Figure6D.pdf")
plot(cpm.mean.pro[, "log2FC"], -log10(as.numeric(cpm.mean.pro[, "adj.P"])), col = cpm.mean.pro[, "color"], 
     pch = 19, cex = 0.3, ylim = c(0,20), xlim = c(-8, 8))
dev.off()
length(which(cpm.mean.pro[, "color"] == "red"))
length(which(cpm.mean.pro[, "color"] == "blue"))


# barplot
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

#dir.create("Figures/")
#pdf("Figures/Figure1B.pdf")
par(mar=c(2,5,1,2))
bp <- barplot(mat, axes = F, beside = T, ylim = c(0,100), names.arg = rep("",3), legend.text = c("-IR","+IR"), space = c(0.2, 0.8))
axis(1, at=c(0,colMeans(bp), colMeans(bp)[3] + 2), labels=F, lwd=3, pos=0)
mtext(text=c("","","","","",""), at=c(0,bp,bp[length(bp)] + 2), line=3, side=1, cex=3)
#axis(2, at = c(0,10,20,30,40,50), cex.axis=2, las=2, lwd=3, labels = T)
axis(2, cex.axis=2, las=2, lwd=3, labels = T)
arrows(bp, mat - mat.sd, bp, mat + mat.sd, angle=90, code=3, length = 0.05, lwd=3)
#dev.off()

