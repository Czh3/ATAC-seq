setwd("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/Rscript")


DNAme = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/DNAme/GSE82234/GSE82234_Matrix_Signal_Intensities.txt", header = T ,row.names = 1, sep = "\t")
GPL13534 = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/DNAme/GSE82234/GPL13534-11288.txt", header = T ,row.names = 1, sep = "\t", quote = "?")


DNAme_beta = data.frame(row.names = rownames(DNAme))
DNAme_beta$donor1_P4 = DNAme$HUVEC_donor1_P4.Methylated.Signal / (DNAme$HUVEC_donor1_P4.Methylated.Signal + DNAme$HUVEC_donor1_P4.Unmethylated.Signal)
DNAme_beta$donor2_P4 = DNAme$HUVEC_donor2_P4.Methylated.Signal / (DNAme$HUVEC_donor2_P4.Methylated.Signal + DNAme$HUVEC_donor2_P4.Unmethylated.Signal)
DNAme_beta$donor3_P4 = DNAme$HUVEC_donor3_P4.Methylated.Signal / (DNAme$HUVEC_donor3_P4.Methylated.Signal + DNAme$HUVEC_donor3_P4.Unmethylated.Signal)
DNAme_beta$donor1_P20 = DNAme$HUVEC_donor1_P20.Methylated.Signal / (DNAme$HUVEC_donor1_P20.Methylated.Signal + DNAme$HUVEC_donor1_P20.Unmethylated.Signal)
DNAme_beta$donor2_P18 = DNAme$HUVEC_donor2_P18.Methylated.Signal / (DNAme$HUVEC_donor2_P18.Methylated.Signal + DNAme$HUVEC_donor2_P18.Unmethylated.Signal)
DNAme_beta$donor3_P13 = DNAme$HUVEC_donor3_P13.Methylated.Signal / (DNAme$HUVEC_donor3_P13.Methylated.Signal + DNAme$HUVEC_donor3_P13.Unmethylated.Signal)

boxplot(DNAme_beta[,1:6])
DNAme_beta = na.omit(DNAme_beta)
DNAme_cor = cor(DNAme_beta[,1:6])
pheatmap(DNAme_cor)
plot(hclust(as.dist(1-cor(DNAme_beta[,1:6]))))

ggplot(DNAme_beta, aes(donor1_P4, donor1_P20)) +
  geom_point( alpha = 0.3, size = 1)+
  #geom_hline(yintercept=0, size = 1, linetype="dashed", color = "gray60") +
  theme_classic(base_size = 20) + #xlab("log10(Distance to TSS)") + ylab("log fold change of ATAC peaks") +
  theme(
    axis.line = element_line(size = 1, colour = "black"),
    axis.text = element_text(colour = "black")) +
  scale_colour_gradientn(colours = brewer.pal(9, "PiYG"))

ggplot(DNAme_beta, aes(donor1_P4, donor2_P4)) +
  geom_point( alpha = 0.3, size = 1)+
  #geom_hline(yintercept=0, size = 1, linetype="dashed", color = "gray60") +
  theme_classic(base_size = 20) + #xlab("log10(Distance to TSS)") + ylab("log fold change of ATAC peaks") +
  theme(
    axis.line = element_line(size = 1, colour = "black"),
    axis.text = element_text(colour = "black")) +
  scale_colour_gradientn(colours = brewer.pal(9, "PiYG"))

## pairs
panel.cor <- function(x, y, digits=2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  test <- cor.test(x,y)
  Signif <- ifelse(round(test$p.value,3)<0.001,"p<0.001",paste("p=",round(test$p.value,3)))  
  text(0.5, 0.25, paste("r=",txt), cex=0.8)
  text(.5, .75, Signif , cex=0.8)
}
panel.smooth<-function (x, y, col = "blue", bg = NA, pch = 18, 
                        cex = 0.8, col.smooth = "red", span = 2/3, iter = 3, ...) 
{
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
          col = col.smooth, ...)
}

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="darkred", ...)
}
png("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/DNAme_paris_scatter.png", res = 300, 1500, 1500)
pairs(DNAme_beta[,1:6], panel= function(...) smoothScatter(..., nrpoints = 0, add = TRUE),
        upper.panel=panel.cor,diag.panel=panel.hist, cex= 0.3)

dev.off() 
  
### DNAme: Distance to TSS
Dis2TSS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/DNAme/GSE82234/methyl.sort.Dis2TSS.bed")
Dis2TSS$Distance = abs(Dis2TSS$V11 - Dis2TSS$V2)
Dis2TSS$young = rowMeans(Dis2TSS[,4:6])
Dis2TSS$old = rowMeans(Dis2TSS[,7:9])
Dis2TSS$FC = log2(Dis2TSS$old/Dis2TSS$young)
#Dis2TSS$Pvalue = apply(Dis2TSS, 1, function(x) wilcox.test(as.numeric(x[4:6]), as.numeric(x[7:9]))$p.value)
#Dis2TSS$Beta =  apply(Dis2TSS, 1, function(x) glm(c(0,0,0,1,1,1) ~ as.numeric(x[4:9]))$coefficients["as.numeric(Dis2TSS[1, 4:9])"])
coefficients["as.numeric(Dis2TSS[1, 4:9])"]

png("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/DNAme_Distance2TSS.png", res = 300, 1500, 1500)
ggplot(na.omit(Dis2TSS), aes(log10(Distance), FC)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept=0, size = 1, linetype="dashed", color = "gray60") +
  theme_classic(base_size = 20) + xlab("log10(Distance to TSS)") + ylab("DNAme: log2(old/young)") +
  theme(
    axis.line = element_line(size = 1, colour = "black"),
    axis.text = element_text(colour = "black")) +
  scale_colour_gradientn(colours = rev(brewer.pal(9, "PiYG")))
dev.off()




# correlation with ATAC
DNAme_atac = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/DNAme/GSE82234/methyl.sort.atac.bed")
DNAme_atac = na.omit(DNAme_atac)
DNAme_atac$V13 = log2(DNAme_atac$V13+0.01)
cor.test(DNAme_atac$V4, DNAme_atac$V13)
par(cex=1.8, lwd=2.5)
png("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/DNAme_ATAC_cor.png", res = 300, 1500, 1500)
smoothScatter(DNAme_atac$V4, DNAme_atac$V13, main="r = -0.68\np-value < 2.2e-16", xlab = "DNA methylation", ylab = "log2( ATAC Signal )")
abline(lm(V13 ~ V4, data = DNAme_atac), col="darkred", lwd=2.5)
dev.off()


DNAme_atac = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/DNAme/GSE82234/methyl.sort.atac1.bed")
DNAme_atac = na.omit(DNAme_atac)
DNAme_atac$Methyl_FC = log2(DNAme_atac$V6 / DNAme_atac$V4)
DNAme_atac$ATAC_FC = log2(DNAme_atac$V17 / DNAme_atac$V13)
DNAme_atac = na.omit(DNAme_atac)
smoothScatter(DNAme_atac$Methyl_FC, DNAme_atac$ATAC_FC, main="r = -0.68\np-value < 2.2e-16", xlab = "DNA methylation", ylab = "log2( ATAC Signal )")

ggplot(DNAme_atac, aes(Methyl_FC, ATAC_FC)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_hline(yintercept=0, size = 1, linetype="dashed", color = "gray60") +
  theme_classic(base_size = 20) + xlab("log10(Distance to TSS)") + ylab("DNAme: log2(old/young)") +
  theme(
    axis.line = element_line(size = 1, colour = "black"),
    axis.text = element_text(colour = "black")) +
  scale_colour_gradientn(colours = rev(brewer.pal(9, "PiYG")))


######
DNAme_beta$chr = GPL13534[rownames(DNAme_beta), "CHR"]
DNAme_beta$position = GPL13534[rownames(DNAme_beta), "MAPINFO"]
#write.table(DNAme_beta, "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/DNAme/GSE82234/methyl.txt", quote = F, sep = "\t", col.names = F )


# IARs DARs
IARs_me = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/DNAme/GSE82234/IARs.methyl.bed")
IARs_me = IARs_me[IARs_me$V4 != ".", ]
IARs_me = data.matrix(IARs_me[,4:9])

pheatmap(IARs_me, scale = "row", color = colorRampPalette(c("blue","black", "red"))(50), show_rownames = F)

DARs_me = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/DNAme/GSE82234/DARs.methyl.bed")
DARs_me = DARs_me[DARs_me$V4 != ".", ]
DARs_me = data.matrix(DARs_me[,4:9])
boxplot(DARs_me)
pheatmap(DARs_me, scale = "row", color = colorRampPalette(c("blue","black", "red"))(50), show_rownames = F)

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/DNAme.I_DARs.pdf", 5, 5)
par(lwd=3)
boxplot(log2(IARs_me[,4]/IARs_me[,1]), log2(IARs_me[,5]/IARs_me[,2]), log2(IARs_me[,6]/IARs_me[,3]),
        log2(DARs_me[,4]/DARs_me[,1]), log2(DARs_me[,5]/DARs_me[,2]), log2(DARs_me[,6]/DARs_me[,3]),
        ylim=c(-2,2), cex = 0.4, col = c("red", "red", "red", "blue", "blue", "blue"))
dev.off()

boxplot(IARs_me)
boxplot(DARs_me)


# FOSL1 motif methylation
withinMotif = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/DNAme/GSE82234/FOSL1.sites.hg19.withinPeaks.methyl.bed")
withinMotif$sites = withinMotif$V8 - withinMotif$V2 + 2
withinMotif[withinMotif$V6 == "-", "sites"] = 10 - withinMotif[withinMotif$V6 == "-", "sites"] + 1
table(withinMotif$sites)

withinMotif.m = melt(withinMotif[,c(10:16)], id = "sites")
ggplot(withinMotif.m, aes(as.factor(sites), value)) + 
  geom_boxplot(width = 0.5, lwd=1, aes(fill=variable)) + #geom_jitter(width = 0.2) +
  theme_classic(base_size = 20) +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"))+
  scale_fill_brewer(palette="Set1") 

withoutMotif = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/DNAme/GSE82234/FOSL1.sites.hg19.withoutPeaks.methyl.bed")
withoutMotif$sites = withoutMotif$V8 - withoutMotif$V2 + 2
withoutMotif[withoutMotif$V6 == "-", "sites"] = 10 - withoutMotif[withoutMotif$V6 == "-", "sites"] +1
table(withoutMotif$sites)

withoutMotif.m = melt(withoutMotif[,c(10:16)], id = "sites")
ggplot(withoutMotif.m, aes(as.factor(sites), value)) + 
  geom_boxplot(width = 0.5, lwd=1, aes(fill=variable)) + #geom_jitter(width = 0.2) +
  theme_classic(base_size = 20) +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"))+
  scale_fill_brewer(palette="Set1") 

withinMotif.m$variable <- "opened"
withoutMotif.m$variable <- "closed"
fols1.motif = rbind(withinMotif.m, withoutMotif.m)
fols1.motif$sites = as.factor(fols1.motif$sites)
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/DNAme_FOSL1_motifs.pdf", 5, 4)
ggplot(fols1.motif, aes(variable, value)) + 
  geom_boxplot(width = 0.5, lwd=1, aes(fill=sites)) + #geom_jitter(width = 0.2) +
  theme_classic(base_size = 20) + ylab("DNA methylation") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "top")+
  scale_fill_brewer(palette="Set1") 
dev.off()
