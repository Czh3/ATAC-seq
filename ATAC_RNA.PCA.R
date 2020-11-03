##### ATAC-seq
# 
library(edgeR)
library(ggplot2)
setwd("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/Rscript")
atac = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/merged_peaks.reads_reform.bed", row.names = 1)
colnames(atac) = c( "midlife.rep2", "old.rep2", "young.rep2", "young.rep1", "midlife.rep1", "old.rep1")
condition = factor(c("midlife", "old", "young", "young", "midlife", "old"), levels = c("young", "midlife", "old"))


#save.image("ATAC.peaksReads.RData")
#load("ATAC.peaksReads.RData")

y <- DGEList(counts=atac,group=condition)
y <- calcNormFactors(y)
# normlized counts
logcpm <- cpm(y, prior.count=1, log=TRUE)
peaks_cpm <- cpm(y, prior.count=1, log=F)

# MA-plot
ggplot(data.frame(peaks_cpm)[1:6000, ], aes(log(old.rep2) + log(old.rep1), log(old.rep2) - log(old.rep1))) +
  geom_point(size=0.5, alpha = 0.3)+
  geom_hline(yintercept = 0) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"))



sds = apply(peaks_cpm, 1, sd)
means = apply(peaks_cpm, 1, mean)
cvs = sds / means
peaks_cpm_selected = peaks_cpm[names(sort(sds, decreasing = T)[1:10000]),]

# pca
condition.info = data.frame("condition" = condition)
row.names(condition.info) = c( "midlife.rep2", "old.rep2", "young.rep2", "young.rep1", "midlife.rep1", "old.rep1")
pca = prcomp(t(peaks_cpm_selected), scale = T)
pca = as.data.frame(pca$x)

pca = cbind(condition, pca)
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/ATAC.PCA.pdf", width = 4.7, height = 3.6)
ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(color = condition), size=4)+
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black")) +
  scale_colour_brewer(palette="Dark2") 
dev.off()

# hclust
d = dist(t(peaks_cpm_selected))
h = hclust(d)
plot(h)


#
design <- model.matrix(~condition)
y <- estimateDisp(y,design)

fit <- glmFit(y, design)

# An ANOVA-like test for any differences
lrt <- glmLRT(fit, coef=3)
sigRegions = topTags(lrt, n=nrow(logcpm))$table
DARs = sigRegions[sigRegions$FDR <= 0.05, ]
DARs= cbind(DARs, data.frame(do.call('rbind', strsplit(rownames(DARs), ":|_")))[2:3])
DARs$len = strtoi(DARs$X3) - strtoi(DARs$X2)
DARs$logCPKM = DARs$logCPM - log(DARs$len / 1000) 

dim(DARs)
dim(DARs[DARs$logCPKM > 0.5,])
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/ATAC.edgeR.FDR0.05.heatmap.pdf", width = 5, height =7)
pheatmap(logcpm[rownames(DARs), ], scale = "row", cluster_rows = T,cluster_cols = T,color = colorRampPalette(c("purple", "black", "yellow"))(50),
         show_rownames = F, show_colnames = T, cutree_rows=2)
dev.off()

write.table(na.omit(DARs), "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/ATAC.DERs.edgeR.FDR0.05.txt", quote = F, sep = "\t")
write.table(na.omit(DARs[DARs$logCPKM > 0.5, -c(6,7)]), "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/ATAC.DERs.edgeR.FDR0.1.txt", quote = F, sep = "\t")
write.table(na.omit(DARs[DARs$logCPKM > 0.5, -c(6,7)]), "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/ATAC.DERs.edgeR.FDR0.05_0.5RPKM.txt", quote = F, sep = "\t")
#write.table(na.omit(sigRegions), "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/ATAC.DERs.edgeR.allPeaks.txt", quote = F, sep = "\t")
write.table(na.omit(sigRegions[(nrow(sigRegions)-999):nrow(sigRegions),]), "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/ATAC.DERs.edgeR.contantPeaks.txt", quote = F, sep = "\t")


########### ATAC peaks annotation ####
ATAC_ann = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/peaks/peak_anno/P5_peaks.bed.ann.variant_function", stringsAsFactors = F)
ATAC_ann = data.frame(table(ATAC_ann$V1))
ATAC_ann[2,2] = ATAC_ann[2,2]+ATAC_ann[5,2]
ATAC_ann[2,2] = ATAC_ann[2,2]+ATAC_ann[8,2]
ATAC_ann[4,2] = ATAC_ann[4,2]+ATAC_ann[6,2]
ATAC_ann[2,2] = ATAC_ann[2,2]+ATAC_ann[7,2]
ATAC_ann = ATAC_ann[-c(5,6,7,8), ]
colnames(ATAC_ann) = c("Region", "Freq")
ATAC_ann$Region = factor(ATAC_ann$Region, levels = c("intergenic", "upstream", "downstream", "exonic", "intronic", "UTR5", "UTR3"))

colors <- c('rgb(211,94,96)', 'rgb(128,133,133)', 'rgb(144,103,167)', 'rgb(171,104,87)', 'rgb(114,147,203)')
m <- list(
  l = 150,
  r = 150,
  b = 200,
  t = 200,
  pad = 4
)
plot_ly(ATAC_ann, labels = ~Region, values = ~Freq, type = 'pie',textfont = list(size = 20),
        textinfo = 'label+percent',
        insidetextfont = list(color = '#FFFFFF'),
        #hoverinfo = 'text',
        marker = list(colors = colors, 
                      line = list(color = '#FFFFFF', width = 1.5)),
        #The 'pull' attribute can also be used to create space between the sectors
        showlegend = T) %>%
  layout( margin = m,
         xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
         yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


# chromstates
ATAC_CS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/allPeaks.ChromState.percent.txt")
ATAC_CS[3,2] = ATAC_CS[3,2] + ATAC_CS[14,2]
ATAC_CS[7,2] = ATAC_CS[7,2] + ATAC_CS[8,2]
ATAC_CS = ATAC_CS[-c(14,8), ]
ATAC_CS = cbind(ATAC_CS, data.frame(do.call('rbind', strsplit(ATAC_CS[,1], "\\d+_", perl = T)))[2])
plot_ly(ATAC_CS, labels = ~X2, values = ~V2, type = 'pie',textfont = list(size = 25),
        textinfo = 'label+percent',
        insidetextfont = list(color = '#FFFFFF'),
        #hoverinfo = 'text',
        #marker = list(colors = colors, 
        #              line = list(color = '#FFFFFF', width = 1.5)),
        #The 'pull' attribute can also be used to create space between the sectors
        showlegend = T) %>%
  layout( margin = m,
          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

ATAC_IARs = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/IARs.ChromState.Percent.txt")
ATAC_IARs[3,2] = ATAC_IARs[3,2] + ATAC_IARs[14,2]
ATAC_IARs[7,2] = ATAC_IARs[7,2] + ATAC_IARs[8,2]
ATAC_IARs = ATAC_IARs[-c(14,8), ]
ATAC_IARs = cbind(ATAC_IARs, data.frame(do.call('rbind', strsplit(ATAC_IARs[,1], "\\d+_", perl = T)))[2])
plot_ly(ATAC_IARs, labels = ~X2, values = ~V2, type = 'pie',textfont = list(size = 20),
        textinfo = 'label+percent',
        insidetextfont = list(color = '#FFFFFF'),
        #hoverinfo = 'text',
        #marker = list(colors = colors, 
        #              line = list(color = '#FFFFFF', width = 1.5)),
        #The 'pull' attribute can also be used to create space between the sectors
        showlegend = T) %>%
  layout( margin = m,
          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))

ATAC_DARs = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/DARs.ChromState.Percent.txt")
ATAC_DARs[3,2] = ATAC_DARs[3,2] + ATAC_DARs[14,2]
ATAC_DARs[7,2] = ATAC_DARs[7,2] + ATAC_DARs[8,2]
ATAC_DARs = ATAC_DARs[-c(14,8), ]
ATAC_DARs = cbind(ATAC_DARs, data.frame(do.call('rbind', strsplit(ATAC_DARs[,1], "\\d+_", perl = T)))[2])
m <- list(
  l = 150,
  r = 200,
  b = 300,
  t = 200,
  pad = 4
)
plot_ly(ATAC_DARs, labels = ~X2, values = ~V2, type = 'pie',textfont = list(size = 25),
        textinfo = 'label+percent',
        insidetextfont = list(color = '#FFFFFF'),
        #hoverinfo = 'text',
        #marker = list(colors = colors, 
        #              line = list(color = '#FFFFFF', width = 1.5)),
        #The 'pull' attribute can also be used to create space between the sectors
        showlegend = T) %>%
  layout( margin = m,
          xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
          yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))


######### ATAC peaks venn plot ####
library(Vennerable)
Vcombo <- Venn(SetNames = c("Young", "Midlife", "old"),
               Weight = c(0, 16233, 1070, 3887, 20245, 32581, 3887, 71872))
plot(Vcombo, show = list(SetLabels = T, Faces = FALSE, FaceText= T))

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/ATAC_peak_venn.pdf", 6, 6)
venn.plot <- draw.triple.venn(
  area1 = 121756,
  area2 = 79847,
  area3 = 128585,
  n12 = 72942,
  n23 = 75759,
  n13 = 104453,
  n123 = 71872,
  category = c("Young", "Midlife", "old"),
  fill = c("darkblue", "darkred", "grey60"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.col = c("darkblue", "darkred", "black")
)
dev.off()


########### distance of ATAC peaks to TSS #####
young_peak = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/peaks/peak_TSS_distence/P5.peak2Tss_dist.bed", stringsAsFactors = F)
qplot(log10(abs(young_peak$V9+1)), binwidth=.05)
mid_peak = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/peaks/peak_TSS_distence/P26.peak2Tss_dist.bed", stringsAsFactors = F)
qplot(log10(abs(mid_peak$V9+1)), binwidth=.05)
old_peak = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/peaks/peak_TSS_distence/P36.peak2Tss_dist.bed", stringsAsFactors = F)
qplot(log10(abs(old_peak$V9+1)), binwidth=.05)
IARs = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/peaks/peak_TSS_distence/IARs.peak2Tss_dist.bed", stringsAsFactors = F)
qplot(log10(abs(IARs$V12+1)), binwidth=.05)
DARs = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/peaks/peak_TSS_distence/DARs.peak2Tss_dist.bed", stringsAsFactors = F)
qplot(log10(abs(DARs$V12+1)), binwidth=.05)

#P4_rep2 = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/peaks/peak_TSS_distence/P4_rep2.peak2Tss_dist.bed", stringsAsFactors = F)
#qplot(log10(abs(P4_rep2$V9+1)), binwidth=.05)

distance_tss = data.frame(rbind(cbind(young_peak$V9, "All Peaks"),
                     cbind(mid_peak$V9, "All Peaks"),
                     cbind(old_peak$V9, "All Peaks"),
                     cbind(IARs$V12, "IARs"),
                     cbind(DARs$V12, "DARs")), stringsAsFactors = F)
colnames(distance_tss) = c("distance", "samples")
distance_tss$distance = as.numeric(distance_tss$distance)
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/ATAC_Distance_peaks_to_TSS.pdf", 7,5)
ggplot(distance_tss, aes(x=log10(abs(distance)+1), fill=samples)) +
  geom_density(alpha=.5, size=1) +
  theme_classic(base_size = 25) + xlab("log10(Distance to TSS)") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"))
dev.off()

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/IARs_DARs_Distance_peaks_to_TSS.pdf", 6,5)
distance_tss$samples = factor(distance_tss$samples, levels = c("All Peaks", "IARs", "DARs"))
ggplot(distance_tss, aes(x=log10(abs(distance)+1), color=samples, fill=samples)) + 
  geom_histogram(aes(y=..count..), alpha=0.5, binwidth=.05, position="identity") +
  facet_grid(samples ~ ., scales = "free_y") +
  geom_vline(xintercept=log10(2000), size = 1, linetype="dashed", color = "gray60") +
  theme_bw(base_size = 25) + xlab("log10(Distance to TSS)") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        legend.position = "none")
dev.off()

wilcox.test(abs(distance_tss[distance_tss$samples == "All Peaks", 1]), abs(distance_tss[distance_tss$samples == "DARs", 1]))
wilcox.test(abs(distance_tss[distance_tss$samples == "All Peaks", 1]), abs(distance_tss[distance_tss$samples == "IARs", 1]))
wilcox.test(abs(distance_tss[distance_tss$samples == "IARs", 1]), abs(distance_tss[distance_tss$samples == "DARs", 1]))

# 
all_peaks = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/peaks/peak_TSS_distence/allPeaks.peak2Tss_dist.bed", stringsAsFactors = F )
colnames(all_peaks) = c("1", "2", "3", "logFC", "5", "6", "7", "FDR", "9", "10", "11", "distance")
qplot(log10(abs(all_peaks$distance)), binwidth=.05)

png("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/ATAC_peaks_distance_to_TSS.scatter.png", width = 450, height = 400)
ggplot(all_peaks, aes(log10(abs(distance)), logFC)) +
  geom_point(aes(color = FDR), alpha = 0.8, size = 1)+
  geom_hline(yintercept=0, size = 1, linetype="dashed", color = "gray60") +
  theme_classic(base_size = 20) + xlab("log10(Distance to TSS)") + ylab("log fold change of ATAC peaks") +
  theme(
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black")) +
  scale_colour_gradientn(colours = brewer.pal(9, "PiYG"))
dev.off()


IARs = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/peaks/peak_TSS_distence/IARs_FDR.peak2Tss_dist.bed", stringsAsFactors = F)
qplot(log10(abs(IARs$V7+1)), binwidth=.05)
DARs = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/20170119/peaks/peak_TSS_distence/DARs_FDR.peak2Tss_dist.bed", stringsAsFactors = F)
qplot(log10(abs(DARs$V7+1)), binwidth=.05)



### ATAC peaks and nearest genes' expression
IARs_gene = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/edgeR.increase.ATAC.peaks.bed.anno", stringsAsFactors = F, sep = "\t", quote = "*", header = T)
DARs_gene = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/edgeR.decrease.ATAC.peaks.bed.anno", stringsAsFactors = F, sep = "\t", quote = "*", header = T)
RNA_expr = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA-seq_Replicate/RNA.cpm.txt", row.names = 1, header = T, stringsAsFactors = F)

RNA_expr_IARs = na.omit(RNA_expr[IARs_gene$Gene.Name, ])

RNA_expr_IARs = melt(RNA_expr_IARs)
RNA_expr_IARs$variable = as.character(RNA_expr_IARs$variable)
RNA_expr_IARs = cbind(RNA_expr_IARs, data.frame(do.call('rbind', strsplit(RNA_expr_IARs$variable, "_"))))
colnames(RNA_expr_IARs) = c("samples", "expr", "state", "rep")
RNA_expr_IARs$samples = factor(RNA_expr_IARs$samples, levels = c("young_rep1","young_rep2", "midlife_rep1","midlife_rep2", "old_rep1","old_rep2"))
RNA_expr_IARs$state = factor(RNA_expr_IARs$state, levels = c("young", "midlife", "old"))
ggplot(RNA_expr_IARs, aes(x=samples, y=log10(expr))) + 
  geom_violin(aes(fill = state)) +
  geom_boxplot(width = 0.1) +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"))+
  scale_fill_brewer(palette="Set1") 

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/IARs_correlate_expr.pdf", 4.5,4)
ggplot(RNA_expr_IARs, aes(x=state, y=log10(expr))) + 
  geom_violin(aes(fill = state), size=1) +
  geom_boxplot(width = 0.1, size=1) +
  geom_hline(yintercept=1.215003, size = 1, linetype="dashed", color = "gray60") +
  theme_classic(base_size = 20) + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none")+
  scale_fill_brewer(palette="Set1") 
dev.off()

summary(RNA_expr_IARs[RNA_expr_IARs$samples %in% c("young_rep1","young_rep2"),])
wilcox.test(RNA_expr_IARs[RNA_expr_IARs$samples %in% c("young_rep1","young_rep2"),2], RNA_expr_IARs[RNA_expr_IARs$samples %in% c("midlife_rep1","midlife_rep2"),2])
wilcox.test(RNA_expr_IARs[RNA_expr_IARs$samples %in% c("midlife_rep1","midlife_rep2"),2], RNA_expr_IARs[RNA_expr_IARs$samples %in% c("old_rep1","old_rep2"),2])
wilcox.test(RNA_expr_IARs[RNA_expr_IARs$samples %in% c("young_rep1","young_rep2"),2], RNA_expr_IARs[RNA_expr_IARs$samples %in% c("old_rep1","old_rep2"),2])


RNA_expr_DARs = na.omit(RNA_expr[DARs_gene$Gene.Name, ])
RNA_expr_DARs = melt(RNA_expr_DARs)
RNA_expr_DARs$variable = as.character(RNA_expr_DARs$variable)
RNA_expr_DARs = cbind(RNA_expr_DARs, data.frame(do.call('rbind', strsplit(RNA_expr_DARs$variable, "_"))))
colnames(RNA_expr_DARs) = c("samples", "expr", "state", "rep")
RNA_expr_DARs$samples = factor(RNA_expr_DARs$samples, levels = c("young_rep1","young_rep2", "midlife_rep1","midlife_rep2", "old_rep1","old_rep2"))
RNA_expr_DARs$state = factor(RNA_expr_DARs$state, levels = c("young", "midlife", "old"))

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/DARs_correlate_expr.pdf", 4.5,4)
ggplot(RNA_expr_DARs, aes(x=state, y=log10(expr))) + 
  geom_violin(aes(fill = state), size=1) +
  geom_boxplot(width = 0.1, size=1) +
  geom_hline(yintercept=1.054958, size = 1, linetype="dashed", color = "gray60") +
  theme_classic(base_size = 20) + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none")+
  scale_fill_brewer(palette="Set1") 
dev.off()

summary(RNA_expr_DARs[RNA_expr_DARs$samples %in% c("young_rep1","young_rep2"),])
wilcox.test(RNA_expr_DARs[RNA_expr_DARs$samples %in% c("young_rep1","young_rep2"),2], RNA_expr_DARs[RNA_expr_DARs$samples %in% c("midlife_rep1","midlife_rep2"),2])
wilcox.test(RNA_expr_DARs[RNA_expr_DARs$samples %in% c("midlife_rep1","midlife_rep2"),2], RNA_expr_DARs[RNA_expr_DARs$samples %in% c("old_rep1","old_rep2"),2])
wilcox.test(RNA_expr_DARs[RNA_expr_DARs$samples %in% c("young_rep1","young_rep2"),2], RNA_expr_DARs[RNA_expr_DARs$samples %in% c("old_rep1","old_rep2"),2])


#bedtools closet gene
IARs_gene = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/edgeR.increase.ATAC.peaks.bed.bedAnn", stringsAsFactors = F, sep = "\t", quote = "*", )
DARs_gene = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/edgeR.decrease.ATAC.peaks.bed.bedAnn", stringsAsFactors = F, sep = "\t", quote = "*", )
expr.up = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA-seq_Replicate/gene.up.txt")
expr.down = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA-seq_Replicate/gene.down.txt")
breaksList = seq(-1.5, 1.5, by = 0.1)

RNA_expr_IARs = na.omit(RNA_expr[intersect(IARs_gene$V15, expr.up$V1), ])
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/IAR_DEGs.pdf", 4, 8)
RNA_expr_IARs = t(scale(t(RNA_expr_IARs), center = T))
RNA_expr_IARs[RNA_expr_IARs >= 1.5] <- 1.5
RNA_expr_IARs[RNA_expr_IARs <= -1.5] <- -1.5
pheatmap(RNA_expr_IARs[,c(3,6,1,5,2,4)], scale = "none", cluster_cols = F, color = colorRampPalette(c("navy", "white", "firebrick3"))(30), border_color = "white")
dev.off()

RNA_expr_DARs = na.omit(RNA_expr[intersect(DARs_gene$V15, expr.down$V1), ])
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/DAR_DEGs.pdf", 4, 8)
RNA_expr_DARs = t(scale(t(RNA_expr_DARs), center = T))
RNA_expr_DARs[RNA_expr_DARs >= 1.5] <- 1.5
RNA_expr_DARs[RNA_expr_DARs <= -1.5] <- -1.5
pheatmap(RNA_expr_DARs[,c(3,6,1,5,2,4)], scale = "none", cluster_cols = F, color = colorRampPalette(c("navy", "white", "firebrick3"))(50), border_color = "white")
dev.off()


RNA_expr_IARs = na.omit(RNA_expr[unique(IARs_gene$V15), ])
pheatmap(RNA_expr_IARs[RNA_expr_IARs$old_rep2 > RNA_expr_IARs$young_rep2 & RNA_expr_IARs$old_rep1 > RNA_expr_IARs$young_rep1, c(3,6,1,5,2,4)], scale = "row", cluster_cols = F)

RNA_expr_IARs = na.omit(RNA_expr[unique(IARs_gene$V15), ])
pheatmap(RNA_expr_IARs, scale = "row")
RNA_expr_IARs = melt(RNA_expr_IARs)
RNA_expr_IARs$variable = as.character(RNA_expr_IARs$variable)
RNA_expr_IARs = cbind(RNA_expr_IARs, data.frame(do.call('rbind', strsplit(RNA_expr_IARs$variable, "_"))))
colnames(RNA_expr_IARs) = c("samples", "expr", "state", "rep")
RNA_expr_IARs$samples = factor(RNA_expr_IARs$samples, levels = c("young_rep1","young_rep2", "midlife_rep1","midlife_rep2", "old_rep1","old_rep2"))
RNA_expr_IARs$state = factor(RNA_expr_IARs$state, levels = c("young", "midlife", "old"))

ggplot(RNA_expr_IARs, aes(x=state, y=log10(expr))) + 
  geom_violin(aes(fill = state), size=1) +
  geom_boxplot(width = 0.1, size=1) +
  geom_hline(yintercept=1.215003, size = 1, linetype="dashed", color = "gray60") +
  theme_classic(base_size = 20) + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none")+
  scale_fill_brewer(palette="Set1") 

# AP1 motif in IARs
IARs_with = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/ENCODE/AP-1/IARs_with_AP1_motif.bed.bedAnn")
IARs_with = na.omit(RNA_expr[unique(IARs_with$V15),])
IARs_with$old = IARs_with$old_rep2 + IARs_with$old_rep1
IARs_with$midlife = IARs_with$midlife_rep2 + IARs_with$midlife_rep1
IARs_with$young = IARs_with$young_rep2 + IARs_with$young_rep1
IARs_with.m = melt(IARs_with[,7:9])
IARs_with.m$variable = factor(IARs_with.m$variable, levels = c("young", "midlife", "old"))
ggplot(IARs_with.m, aes(x=variable, y=log10(value))) + 
  geom_violin(aes(fill = variable), size=1) +
  geom_boxplot(width = 0.1, size=1) +
  #geom_hline(yintercept=1.215003, size = 1, linetype="dashed", color = "gray60") +
  theme_classic(base_size = 20) + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none")+
  scale_fill_brewer(palette="Set1") 



IARs_without = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/ENCODE/AP-1/IARs_without_AP1_motif.bed.bedAnn")
IARs_without = na.omit(RNA_expr[unique(IARs_without$V15),])
IARs_without$old = IARs_without$old_rep2 + IARs_without$old_rep1
IARs_without$midlife = IARs_without$midlife_rep2 + IARs_without$midlife_rep1
IARs_without$young = IARs_without$young_rep2 + IARs_without$young_rep1
IARs_without.m = melt(IARs_without[,7:9])
IARs_without.m$variable = factor(IARs_without.m$variable, levels = c("young", "midlife", "old"))
ggplot(IARs_without.m, aes(x=variable, y=log10(value))) + 
  geom_violin(aes(fill = variable), size=1) +
  geom_boxplot(width = 0.1, size=1) +
  #geom_hline(yintercept=1.215003, size = 1, linetype="dashed", color = "gray60") +
  theme_classic(base_size = 20) + xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none")+
  scale_fill_brewer(palette="Set1") 


# IARs in heterochrom/lo
IARs_Hetero_gene = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/IARs.Heterochrom.bed.bedAnn", stringsAsFactors = F, sep = "\t", quote = "*", )
  
expr = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA-seq_Replicate/RNA.cpm.txt")
IARs_Hetero_gene_expr = na.omit(unique(expr[IARs_Hetero_gene$V15, ]))
pheatmap(IARs_Hetero_gene_expr[rowMeans(IARs_Hetero_gene_expr) >= 1, c(3,6,1,5,2,4)], scale = "row", cluster_cols = F)
boxplot(log2(IARs_Hetero_gene_expr[, c(3,6,1,5,2,4)])+1)

# conservation hg19.100way.phastCons.bw
Random = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/Cons/random.tab")
Random = data.frame(Cons = as.numeric(Random$V5),
                  Sample = "Random")
AllPeaks = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/Cons/merged_peaks.tab")
AllPeaks = data.frame(Cons = as.numeric(AllPeaks$V5),
                    Sample = "AllPeaks")
IARs_Het = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/Cons/IARs.Heterochrom.tab")
IARs_Het = data.frame(Cons = as.numeric(IARs_Het$V5),
                    Sample = "IARs_Het")

df = rbind(Random, AllPeaks, IARs_Het)
df$Sample = factor(df$Sample, levels = c("Random", "AllPeaks", "IARs_Het"))
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/Conservation_peaks.pdf", 6, 4)
ggplot(df, aes(x=Cons)) +
  geom_density(aes(group=Sample, colour=Sample), size = 1) +
  theme_classic(base_size = 20) +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black")) +
  xlab("100way.phastCons") + #xlim(0,0.6) +
  scale_colour_brewer(palette="Set1")
dev.off()
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/Conservation_peaks_zoomin.pdf", 10, 4)
ggplot(df, aes(x=Cons)) +
  geom_density(aes(group=Sample, colour=Sample), size = 0.8) +
  theme_classic(base_size = 20) +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none") +
  xlab("100way.phastCons") + xlim(0,0.20) + ylim(2.5,10) +
  scale_colour_brewer(palette="Set1")
dev.off()

# IARs ENCODE:DNase 
IAR_DNase = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/IARs.Heterochrom.3.bed.412")
DNase_name = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/ENCODE/DNase/names.txt")
IAR_DNase = IAR_DNase[,-c(1:3)]
colnames(IAR_DNase) = 1:ncol(IAR_DNase)
IAR_DNase = apply(IAR_DNase, 2, function(x) ifelse(x>=1, 1 ,0))
DNase_name = DNase_name[,1:2]
colnames(DNase_name) = c("Tissue", "Age")
Tissue = data.frame(DNase_name[DNase_name$Age == "fetal", 1])
rownames(Tissue) = rownames(DNase_name[DNase_name$Age == "fetal", ])
colnames(Tissue) = c("Tissue")

png("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/IARs_Hetero_Roadmap.DNase.png", 800, 800)
pheatmap(t(IAR_DNase[,rownames(Tissue)]), cluster_cols=T, annotation_row = Tissue,
         cluster_rows=F,  show_rownames=F, show_colnames=F,fontsize=12, clustering_method="ward.D", #clustering_distance_rows = "binary",clustering_distance_cols = "binary",
         color = colorRampPalette(c("gray", "firebrick3"))(2))
dev.off()

## ENCODE chrom state
young_CS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks/P4_peaks.bed.chromStat")
midlife_CS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks/P14_peaks.bed.chromStat")
ols_CS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks/P26_peaks.bed.chromStat")

chromState = cbind(young_CS, midlife_CS, ols_CS)
chromState = chromState[,c(1,2,4,6)]
colnames(chromState) = c("chromState", "Young", "Midlife", "Old")
chromState.m = melt(chromState)
chromState.m$chromState = factor(chromState.m$chromState, 
                                 levels = c("15_Repetitive/CNV","14_Repetitive/CNV", "13_Heterochrom/lo", "12_Repressed",
                                            "11_Weak_Txn", "10_Txn_Elongation", "9_Txn_Transition", "8_Insulator", "7_Weak_Enhancer",
                                            "6_Weak_Enhancer", "5_Strong_Enhancer", "4_Strong_Enhancer", "3_Poised_Promoter", "2_Weak_Promoter", "1_Active_Promoter"))

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/ATAC_peaks_ChromState.pdf")
ggplot(data=chromState.m, aes(x=chromState , y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  coord_flip() +
  geom_hline(yintercept=0, size = 1, color = "black") +
  theme_classic(base_size = 20) + ylab("log ratio of observed to random") + xlab("Chrom State") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.line.y = element_line(size = 0, colour = "black"),
        legend.position = "top")
dev.off()


IARs_CS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/IARs.ChromState.Percent.txt")
DARs_CS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/DARs.ChromState.Percent.txt")
AllPeaks_CS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/AllPeaks.ChromState.Percent.txt")

chromState = cbind(IARs_CS, DARs_CS, AllPeaks_CS)
chromState = chromState[,c(1,2,4,6)]
colnames(chromState) = c("chromState", "IARs", "DARs", "AllPeaks")
chromState.m = melt(chromState)
chromState.m$chromState = factor(chromState.m$chromState, 
                                 levels = rev(c("15_Repetitive/CNV","14_Repetitive/CNV", "13_Heterochrom/lo", "12_Repressed",
                                            "11_Weak_Txn", "10_Txn_Elongation", "9_Txn_Transition", "8_Insulator", "7_Weak_Enhancer",
                                            "6_Weak_Enhancer", "5_Strong_Enhancer", "4_Strong_Enhancer", "3_Poised_Promoter", "2_Weak_Promoter", "1_Active_Promoter")))
chromState.m$variable = factor(chromState.m$variable, levels = c("AllPeaks", "IARs", "DARs"))

ggplot(data=chromState.m, aes(x=chromState , y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #coord_flip() +
  #geom_hline(yintercept=0, size = 1, color = "black") + 
  theme_classic(base_size = 20) + ylab("log ratio of observed to espected") + xlab("Chrom State") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.text = element_text(colour = "black"),
        #axis.line.y = element_line(size = 0, colour = "black"),
        legend.position = "top")


IARs_CS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/FDR0.05_RPKM0.5/IARs.ChromState.txt")
DARs_CS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/FDR0.05_RPKM0.5/DARs.ChromState.txt")

chromState = cbind(IARs_CS, DARs_CS)
rownames(chromState) = chromState$V1
chromState = chromState[,c(2,4)]
colnames(chromState) = c("IARs", "DARs")
chromState = chromState[c("15_Repetitive/CNV","14_Repetitive/CNV", "13_Heterochrom/lo", "12_Repressed",
                          "11_Weak_Txn", "10_Txn_Elongation", "9_Txn_Transition", "8_Insulator", "7_Weak_Enhancer",
                          "6_Weak_Enhancer", "5_Strong_Enhancer", "4_Strong_Enhancer", "3_Poised_Promoter", "2_Weak_Promoter", "1_Active_Promoter"), ]

pheatmap(chromState, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("white", "black"))(30))


chromState.m = melt(data.matrix(chromState))
chromState.m$chromState = factor(chromState.m$X1, 
                                 levels = rev(c("15_Repetitive/CNV","14_Repetitive/CNV", "13_Heterochrom/lo", "12_Repressed",
                                                "11_Weak_Txn", "10_Txn_Elongation", "9_Txn_Transition", "8_Insulator", "7_Weak_Enhancer",
                                                "6_Weak_Enhancer", "5_Strong_Enhancer", "4_Strong_Enhancer", "3_Poised_Promoter", "2_Weak_Promoter", "1_Active_Promoter")))
chromState.m$variable = factor(chromState.m$X2, levels = c("AllPeaks", "IARs", "DARs"))

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/I_DARs_peaks_ChromState_obs2exp.pdf", 7, 6)
ggplot(data=chromState.m, aes(x=chromState , y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #coord_flip() +
  geom_hline(yintercept=0, size = 1, color = "gray") + 
  theme_classic(base_size = 20) + ylab("log ratio of observed to expected") + xlab("Chrom State") +
  #scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.text = element_text(colour = "black"),
        #axis.line.y = element_line(size = 0, colour = "black"),
        legend.position = "top")
dev.off()

### repeat mask
IARs_RM = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/repeatMask/IARs.enrich.txt")
DARs_RM = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/repeatMask/DARs.enrich.txt")

RepeatMask = cbind(IARs_RM, DARs_RM)
rownames(RepeatMask) = RepeatMask$V1
RepeatMask = RepeatMask[,c(2,4)]
colnames(RepeatMask) = c("IARs", "DARs")
RepeatMask = RepeatMask[-c(3,7,8,12,13,14,15,20), ]

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/I_DARs_peaks_RepeatMask_family_obs2exp.pdf", 3, 6)
brk = c(seq(-6.2,0,length.out = 10), seq(0.2,2,length.out = 9))
pheatmap(RepeatMask, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("blue", "black", "yellow"))(length(brk)), breaks = brk, border_color = NA, fontsize = 15)
dev.off()

#sub-family
IARs_RM = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/repeatMask/IARs.enrich_sub.txt")
DARs_RM = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/repeatMask/DARs.enrich_sub.txt")
RepeatMask = cbind(IARs_RM, DARs_RM)
rownames(RepeatMask) = RepeatMask$V1
RepeatMask = RepeatMask[,c(2,4)]
colnames(RepeatMask) = c("IARs", "DARs")
RepeatMask = RepeatMask[c(1,7,8,11,26,30,32,35,45,51,53), ]

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/I_DARs_peaks_RepeatMask_sub-family_obs2exp.pdf", 3, 6)
brk = seq(-2, 2, 0.1)
pheatmap(RepeatMask, cluster_rows = F, cluster_cols = F, color = colorRampPalette(c("blue", "black", "yellow"))(length(brk)), breaks = brk, border_color = NA, fontsize = 15)
dev.off()

#### histone modification
hm = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/ENCODE/I_DARs.histone.tab", sep="\t", skip=2, stringsAsFactors = F)
colnames(hm) = c("modification", "region", 1:600)
hm.m = melt(hm)
hm.m[hm.m$region == "merged_peaks.bed", "region"] = "All peaks"
hm.m[hm.m$region == "edgeR.increase.ATAC.peaks.bed", "region"] = "IARs"
hm.m[hm.m$region == "edgeR.decrease.ATAC.peaks.bed", "region"] = "DARs"
hm.m$region = factor(hm.m$region, levels = c("All peaks", "IARs", "DARs"))
hm.m$variable = as.numeric(hm.m$variable)
hm.m$modification = do.call('rbind', strsplit(hm.m$modification, "[.]"))[,1]

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/Peak_encode_histone_modification.pdf", 5,16)
p <- ggplot(hm.m, aes(x=variable, y=value, col=region)) + geom_line(size=2) +
  #geom_vline(xintercept = 2000, colour="#E69F00", size=1.2,alpha=0.8) +
  facet_grid(modification ~ ., scales = "free") +
  theme_bw(base_size=25)+
  theme(axis.line = element_line(size = 1, colour = "black"),
        legend.position = "top") +
  labs (y='Chip-seq signal') +
  scale_color_brewer(palette="Set1") 
p + scale_x_continuous(name="", limits=c(0, 600),breaks=c(0, 300, 600), labels=c('-3k','peaks center','3k'))
dev.off()

hm.m1 = hm.m[hm.m$modification %in% c("H3K27ac", "H3K4me3"), ]
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/Peak_encode_H3K27ac_H3K4me3.pdf", 5, 6)
p <- ggplot(hm.m1, aes(x=variable, y=value, col=region)) + geom_line(size=2) +
  #geom_vline(xintercept = 2000, colour="#E69F00", size=1.2,alpha=0.8) +
  facet_grid(modification ~ ., scales = "free") +
  theme_bw(base_size=25)+
  theme(axis.line = element_line(size = 1, colour = "black"),
        legend.position = "top") +
  labs (y='Chip-seq signal') +
  scale_color_brewer(palette="Set1") 
p + scale_x_continuous(name="", limits=c(0, 600),breaks=c(0, 300, 600), labels=c('-3k','peaks center','3k'))
dev.off()


#### motif enrichment
IARs_MEME = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/motif/IARs.MEME/centrimo.txt", quote = "#", stringsAsFactors = F)
IARs_MEME = IARs_MEME[,c(2,6)]
IARs_MEME$V2 = do.call('rbind', strsplit(IARs_MEME$V2, "_"))[,1]

DARs_MEME = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/motif/DARs.MEME/centrimo.txt", quote = "#", stringsAsFactors = F)
DARs_MEME = DARs_MEME[,c(2,6)]
DARs_MEME$V2 = do.call('rbind', strsplit(DARs_MEME$V2, "_"))[,1]

IARs_Homer = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/motif/IARsdistal_motif_unmark/knownResults.txt",sep="\t", stringsAsFactors = F, skip=1)
IARs_Homer = IARs_Homer[,c(1,3)]
IARs_Homer$V1 = do.call('rbind', strsplit(IARs_Homer$V1, "[(]"))[,1]
IARs_Homer$V1 = toupper(IARs_Homer$V1)
IARs_Homer[IARs_Homer$V1 == "FRA1", "V1"] = "FOSL1"

DARs_Homer = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/motif/DARsdistal_motif_unmark/knownResults.txt",sep="\t", stringsAsFactors = F, skip=1)
DARs_Homer = DARs_Homer[,c(1,3)]
DARs_Homer$V1 = do.call('rbind', strsplit(DARs_Homer$V1, "[(]"))[,1]
DARs_Homer$V1 = toupper(IARs_Homer$V1)
DARs_Homer[DARs_Homer$V1 == "FRA1", "V1"] = "FOSL1"

motifs = intersect(intersect(IARs_MEME$V2, DARs_MEME$V2), intersect(IARs_Homer$V1, DARs_Homer$V1))
IARs_Homer = IARs_Homer[IARs_Homer$V1 %in% motifs, ]
DARs_Homer = DARs_Homer[DARs_Homer$V1 %in% motifs, ]
IARs_MEME = IARs_MEME[IARs_MEME$V2 %in% motifs,]
DARs_MEME = DARs_MEME[DARs_MEME$V2 %in% motifs,]
IARs_Homer = IARs_Homer[!duplicated(IARs_Homer[,1]), ]
DARs_Homer = DARs_Homer[!duplicated(DARs_Homer[,1]), ]
motifs_enrich = cbind(IARs_MEME[order(IARs_MEME$V2), ],
                      DARs_MEME[order(DARs_MEME$V2), ],
                      IARs_Homer[order(IARs_Homer$V1), ],
                      DARs_Homer[order(DARs_Homer$V1), ])

motifs_enrich$Homer = -log10(motifs_enrich[,6]) + log10(motifs_enrich[,8])
motifs_enrich$MEME = -log10(motifs_enrich[,2]) + log10(motifs_enrich[,4])
motifs_enrich$mean = rowMeans(motifs_enrich[, 9:10])
motifs_enrich$name = motifs_enrich$V1
motifs_enrich[abs(motifs_enrich$mean) <= 3, "name"] = ""

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/motif_homer_MEME.pdf", 5, 5)
ggplot(motifs_enrich, aes(MEME, Homer, name = V1  )) +
  geom_vline(xintercept = 0, colour="gray60", size=1,alpha=0.8) + 
  geom_hline(yintercept = 0, colour="gray60", size=1,alpha=0.8) +
  geom_point(aes(size = abs(mean), col = log2(abs(mean))))+
  geom_text(aes(label=name), size=4, hjust=-0.2, vjust=0.5) +
  theme_classic(base_size = 25) + xlim(-10,12)+
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"),
        legend.position = "none") +
  scale_colour_distiller(palette="RdGy")
dev.off()

### add constant
IARs_MEME = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/motif/IARs.MEME/centrimo.txt", quote = "#", stringsAsFactors = F)
IARs_MEME = IARs_MEME[,c(2,6)]
IARs_MEME$V2 = do.call('rbind', strsplit(IARs_MEME$V2, "_"))[,1]

DARs_MEME = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/motif/DARs.MEME/centrimo.txt", quote = "#", stringsAsFactors = F)
DARs_MEME = DARs_MEME[,c(2,6)]
DARs_MEME$V2 = do.call('rbind', strsplit(DARs_MEME$V2, "_"))[,1]

Constant_MEME = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/motif/contantPeaks.MEME/centrimo.txt", quote = "#", stringsAsFactors = F)
Constant_MEME = Constant_MEME[,c(2,6)]
Constant_MEME$V2 = do.call('rbind', strsplit(Constant_MEME$V2, "_"))[,1]

MEME_motif = cbind(Constant_MEME, IARs_MEME, DARs_MEME)
MEME_motif$m = rowMeans(MEME_motif[,c(4,6)])
MEME_motif = MEME_motif[order(MEME_motif$m), ]
MEME_motif = MEME_motif[1:20, ]
rownames(MEME_motif) = MEME_motif[,1]
MEME_motif = MEME_motif[,c(1, 2,4,6)]
colnames(MEME_motif) = c("TFs", "Constant", "IARs", "DARs")
MEME_motif.m = melt(MEME_motif)
MEME_motif.m$value = -log10(MEME_motif.m$value)

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/motif_enrich_MEME.pdf", 3, 5)
ggplot(MEME_motif.m, aes(x = factor(variable), y = TFs)) +
  geom_point(aes(fill= value , size = value), shape = 21) +
  scale_fill_gradient(low="blue", high = "gold") +
  theme(text = element_text( colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(size = 13,angle = 45, vjust = 1, hjust = 1),
        axis.ticks = element_line(size = 0),
        #legend.position = "bottom",
        panel.background = element_blank() ) +
  xlab("") + ylab("") + labs(size="-log10(p value)", fill="-log10(p value)")
dev.off()

IARs_Homer = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/motif/IARsdistal_motif_unmark/knownResults.txt",sep="\t", stringsAsFactors = F, skip=1)
IARs_Homer = IARs_Homer[,c(1,3)]
IARs_Homer$V1 = do.call('rbind', strsplit(IARs_Homer$V1, "[(]"))[,1]

DARs_Homer = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/motif/DARsdistal_motif_unmark/knownResults.txt",sep="\t", stringsAsFactors = F, skip=1)
DARs_Homer = DARs_Homer[,c(1,3)]
DARs_Homer$V1 = do.call('rbind', strsplit(DARs_Homer$V1, "[(]"))[,1]

Constant_Homer = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/motif/contantPeaksdistal_motif_unmark/knownResults.txt",sep="\t", stringsAsFactors = F, skip=1)
Constant_Homer = Constant_Homer[,c(1,3)]
Constant_Homer$V1 = do.call('rbind', strsplit(Constant_Homer$V1, "[(]"))[,1]

Homer_motif = cbind(Constant_Homer[order(Constant_Homer$V1),], IARs_Homer[order(IARs_Homer$V1),], DARs_Homer[order(DARs_Homer$V1),])
Homer_motif$m = rowMeans(Homer_motif[,c(4,6)])
Homer_motif = Homer_motif[order(Homer_motif$m), ]
Homer_motif = Homer_motif[1:20, ]
rownames(Homer_motif) = Homer_motif[,1]
Homer_motif = Homer_motif[,c(1, 2,4,6)]
colnames(Homer_motif) = c("TFs", "Constant", "IARs", "DARs")
Homer_motif.m = melt(Homer_motif)
Homer_motif.m$value = -log10(Homer_motif.m$value)

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/motif_enrich_homer.pdf", 3.6, 5)
ggplot(Homer_motif.m, aes(x = factor(variable), y = TFs)) +
  geom_point(aes(fill= value , size = value), shape = 21) +
  scale_fill_gradient( low="green", high = "purple") +
  theme(text = element_text( colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        axis.text.x = element_text(size = 13,angle = 45, vjust = 1, hjust = 1),
        axis.ticks = element_line(size = 0),
        #legend.position = "bottom",
        panel.background = element_blank() ) +
  xlab("") + ylab("") + labs(size="-log10(p value)", fill="-log10(p value)")
dev.off()

#
# FOSL1 motif
FOSL1_motif = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/ENCODE/I_DARs.FOSL1.ATAC.tab",sep="\t", stringsAsFactors = F, skip=2)

FOSL1_motif[1, 3:402] = colMeans(FOSL1_motif[c(1,7), 3:402])
FOSL1_motif[2, 3:402] = colMeans(FOSL1_motif[c(2,8), 3:402])
FOSL1_motif[3, 3:402] = colMeans(FOSL1_motif[c(3,9), 3:402])
FOSL1_motif[4, 3:402] = colMeans(FOSL1_motif[c(4,10), 3:402])
FOSL1_motif[5, 3:402] = colMeans(FOSL1_motif[c(5,11), 3:402])
FOSL1_motif[6, 3:402] = colMeans(FOSL1_motif[c(6,12), 3:402])
FOSL1_motif = FOSL1_motif[1:6,]

colnames(FOSL1_motif) = c("sample", "region", 1:400)
FOSL1_motif.m = melt(FOSL1_motif)
FOSL1_motif.m[FOSL1_motif.m$region == "FOSL1.sites.IARs.bed", "region"] = "IARs"
FOSL1_motif.m[FOSL1_motif.m$region == "FOSL1.sites.DARs.bed", "region"] = "DARs"
FOSL1_motif.m$region = factor(FOSL1_motif.m$region, levels = c("IARs", "DARs"))
FOSL1_motif.m$variable = as.numeric(FOSL1_motif.m$variable)
FOSL1_motif.m$sample = do.call('rbind', strsplit(FOSL1_motif.m$sample, "_"))[,1]
FOSL1_motif.m$sample = factor(FOSL1_motif.m$sample, levels = c("young", "midlife", "old"))

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/I_DARs_FOSL1.motif.pdf", 5,16)
p <- ggplot(FOSL1_motif.m, aes(x=variable, y=value, col=sample)) + geom_line( size=2) +
  #geom_vline(xintercept = 2000, colour="#E69F00", size=1.2,alpha=0.8) +
  facet_grid(region ~ ., scales = "free") +
  theme_bw(base_size=25)+
  theme(axis.line = element_line(size = 1, colour = "black"),
        legend.position = "top") +
  labs (y='Chip-seq signal') +
  scale_color_brewer(palette="Set1") 
p + scale_x_continuous(name="", limits=c(0, 400),breaks=c(0, 200, 400), labels=c('-2k','peaks center','2k'))
dev.off()

# IARs
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/IARs_FOSL1.motif.pdf", 5,6)
p <- ggplot(FOSL1_motif.m[FOSL1_motif.m$region == "IARs", ], aes(x=variable, y=value, col=sample)) + geom_line( size=2) +
  #geom_vline(xintercept = 2000, colour="#E69F00", size=1.2,alpha=0.8) +
  theme_bw(base_size=25)+
  theme(panel.border = element_rect(size = 1, colour = "black", fill = NA), 
        legend.position = "top") +
  labs (y='RPKM (ATAC)') +
  scale_color_brewer(palette="Set1") 
p + scale_x_continuous(name="", limits=c(0, 400),breaks=c(0, 200, 400), labels=c('-2k','peaks center','2k'))
dev.off()
# DARs
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/DARs_FOSL1.motif.pdf", 5,6)
p <- ggplot(FOSL1_motif.m[FOSL1_motif.m$region == "DARs", ], aes(x=variable, y=value, col=sample)) + geom_line( size=2) +
  #geom_vline(xintercept = 2000, colour="#E69F00", size=1.2,alpha=0.8) +
  theme_bw(base_size=25)+
  theme(panel.border = element_rect(size = 1, colour = "black", fill = NA), 
        legend.position = "top") +
  labs (y='RPKM (ATAC)') +
  scale_color_brewer(palette="Set1") 
p + scale_x_continuous(name="", limits=c(0, 400),breaks=c(0, 200, 400), labels=c('-2k','peaks center','2k'))
dev.off()



####### GC content
IARs_GC = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/GC_content/IARs.GC")
DARs_GC = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/GC_content/DARs.GC")
Constant_GC = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/GC_content/Constant.GC")

boxplot(IARs_GC$V10, DARs_GC$V10, Constant_GC$V5)




### insert size
B_pro = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/OIS/ATAC_new/mapping/B-proliferation.sort.rmdup.bam.insert_size_metrics.txt", skip = 11, row.names = 1)
B_ras = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/OIS/ATAC_new/mapping/B-RAS.sort.rmdup.bam.insert_size_metrics.txt", skip = 11, row.names = 1)
D_pro = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/OIS/ATAC_new/mapping/D-proliferation.sort.rmdup.bam.insert_size_metrics.txt", skip = 11, row.names = 1)
D_ras = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/OIS/ATAC_new/mapping/D-RAS.sort.rmdup.bam.insert_size_metrics.txt", skip = 11, row.names = 1)
E_pro = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/OIS/ATAC_new/mapping/E-proliferation.sort.rmdup.bam.insert_size_metrics.txt", skip = 11, row.names = 1)
E_ras = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/OIS/ATAC_new/mapping/E-RAS.sort.rmdup.bam.insert_size_metrics.txt", skip = 11, row.names = 1)

insert_size = data.frame(cbind(B_pro[30:500,1], B_ras[30:500,1], D_pro[30:500,1], D_ras[30:500,1], E_pro[30:500,1], E_ras[30:500,1]))
insert_size[,1] = insert_size[,1]/colSums(insert_size)[1]
insert_size[,2] = insert_size[,2]/colSums(insert_size)[2]
insert_size[,3] = insert_size[,3]/colSums(insert_size)[3]
insert_size[,4] = insert_size[,4]/colSums(insert_size)[4]
insert_size[,5] = insert_size[,5]/colSums(insert_size)[5]
insert_size[,6] = insert_size[,6]/colSums(insert_size)[6]

pheatmap(t(insert_size), show_colnames = T, cluster_col= F, color = colorRampPalette(c("white", "black"))(10))


### chip-seq

IARs_tf = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/motif/IARs.bed.341")
IARs_tf = IARs_tf[, 4:344]
name = read.table("/lustre/user/liclab/zhangc/Taolab/Czh/ENCODE_data/ChIP-seq/K562/download/expr/names.txt")
IARs_tf = apply(IARs_tf, 2, function(x) ifelse(x>=1, 1 ,0))
colnames(IARs_tf) = name$V1

pheatmap(t(IARs_tf), cluster_cols=T, cluster_rows=T, show_rownames=T, show_colnames=F, fontsize=6,
         color = colorRampPalette(c("gray", "firebrick3"))(50))

DARs_tf = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/motif/DARs.bed.341")
DARs_tf = DARs_tf[, 4:344]
DARs_tf = apply(DARs_tf, 2, function(x) ifelse(x>=1, 1 ,0))
colnames(DARs_tf) = name$V1

pheatmap(t(DARs_tf), cluster_cols=T, cluster_rows=T, show_rownames=T, show_colnames=F, fontsize=6,
         color = colorRampPalette(c("gray", "firebrick3"))(50))

