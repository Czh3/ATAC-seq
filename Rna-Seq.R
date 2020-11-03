# RNA-seq for HEVUC senescence project
library(DESeq)
library(edgeR)
library(gplots)
library(reshape)
rep1 = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA/pipe_out/cufflinks/cuffnorm/genes.count_table",  header = T, row.names = 1)
rep2 = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA_rep2/pipe_out/cufflinks/cuffnorm/genes.count_table", header = T, row.names = 1)
rep1 = rep1[, -3]
expr = merge(rep1, rep2, by = "row.names")
rownames(expr) = expr$Row.names
expr = expr[,-1]
expr = floor(expr)
colnames(expr) = c("midlife_rep1", "old_rep1", "young_rep1", "old_rep2", "midlife_rep2", "young_rep2")

condition = c("midlife", "old", "young", "old", "midlife", "young")

y <- DGEList(counts=expr, group=condition)
y <- calcNormFactors(y)
#filter
y <- y[rowSums(cpm(y)>1) >= 2, , keep.lib.sizes=FALSE]

y2 = y
y2$samples$group <- relevel(y2$samples$group, ref="young")

# normlized counts
logcpm <- cpm(y, prior.count=1, log=TRUE)
gene_cpm <- cpm(y, prior.count=1, log=F)

# PCA
sds = apply(gene_cpm[rowSums(gene_cpm) >= 6, ], 1, sd)
means = apply(gene_cpm[rowSums(gene_cpm) >= 6, ], 1, mean)
cvs = sds / means
gene_cpm_selected = gene_cpm[names(sort(cvs, decreasing = T)[1:5000]),]

pca = prcomp(t(gene_cpm), scale = T)
pca = as.data.frame(pca$x)

pca = cbind(condition, pca)
pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/RNA.PCA.pdf", width = 4.7, height = 3.6)
ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(color = condition), size=4)+
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black")) +
  scale_colour_brewer(palette="Dark2") 
dev.off()
# hclust
d = dist(t(logcpm))
h = hclust(d)
plot(h)

save.image(file="HEVUC.RNA.RData")
load("HEVUC.RNA.RData")
write.table(na.omit(gene_cpm), "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA-seq_Replicate/RNA.cpm.txt", quote = F, sep = "\t")


design <- model.matrix(~condition)
y <- estimateDisp(y, design)
fit <- glmFit(y, design)

design2 <- model.matrix(~group, data=y2$samples)
y2 <- estimateDisp(y2, design2)
fit2 <- glmFit(y2, design2)


# An ANOVA-like test for any differences
lrt <- glmLRT(fit, coef=2:3)
sigRegions = topTags(lrt, n=nrow(expr))$table
DEGs = sigRegions[sigRegions$FDR <= 0.05, ]

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/RNA.DEGs.edgeR.heatmap.pdf", width = 5, height =7)
pheatmap(logcpm[rownames(DEGs), ], scale = "row", cluster_rows = T,cluster_cols = T,color = colorRampPalette(c("black", "white", "firebrick3"))(50),
         show_rownames = F, show_colnames = T, cutree_rows=2)
dev.off()

write.table(na.omit(DEGs), "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA-seq_Replicate/RNA.DEGs.edgeR.FDR0.05.txt", quote = F, sep = "\t")
write.table(na.omit(sigRegions), "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA-seq_Replicate/RNA.DEGs.edgeR.all.txt", quote = F, sep = "\t")

# for GSEA rank file
lrt <- glmLRT(fit, coef=3)
sigRegions = topTags(lrt, n=nrow(expr))$table
sigRegions$name = rownames(sigRegions)
sigRegions$rank = -log10(sigRegions$PValue) * sign(sigRegions$logFC)

write.table(na.omit(sigRegions[, c(6,7)]), "/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA-seq_Replicate/RNA.DEGs.edgeR.all.rnk", quote = F, sep = "\t", row.names = F)

# DEGs scatter plot
expr_DEGs = cbind(gene_cpm[rownames(sigRegions), ], sigRegions)
expr_DEGs$FC = log2((expr_DEGs$young_rep1 + expr_DEGs$young_rep2)/(expr_DEGs$midlife_rep1 + expr_DEGs$midlife_rep2))
expr_DEGs$Mean = rowMeans(expr_DEGs[, c(2,3,4,6)])
p = ggplot(expr_DEGs, aes(log(Mean), FC)) +
  geom_point(aes(color = PValue, a=rownames(expr_DEGs)), alpha = 0.8, size = 1) + 
  scale_colour_gradientn(colours = rev(brewer.pal(9, "YlOrRd"))) +
  theme_classic(base_size = 25) +
  xlab("log(mean of expression)") + ylab("log2( Young / Old )") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black"))

png("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/RNA_Young-Old.scatter.png")
p
dev.off()

expr_DEGs$lab = rownames(expr_DEGs)
select.gene = c("SULF1", "GRIK5","ALDH1A1", "IRF6", "BGN", "MKI67", "CDH11",
                "CDKN2A", "CD44","PLAT", "FGF5","KCNMA1", "IGFBP6", "TGFBI")

ggplot(expr_DEGs, aes(log(Mean), FC)) +
  geom_point(aes(color = PValue), alpha = 0.8, size = 1) + 
  scale_colour_gradientn(colours = rev(brewer.pal(9, "YlOrRd"))) +
  theme_classic(base_size = 20) +
  xlab("log(mean of expression)") + ylab("log2( Young / Old )") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black")) +
  geom_text_repel(data = expr_DEGs[expr_DEGs$lab %in% select.gene, ],
                  aes(label = lab),
                  size=5,
                  box.padding = unit(0.35, "lines"),
                  point.padding = unit(0.3, "lines"))

ggsave("../figures/RNA_Young-Old.scatter.pdf", device = "pdf", width = 8, height = 6)


library(plotly)
#RNA_Young-Old.scatter_ggplotly.html
ggplotly(p)

ph = pheatmap(logcpm[rownames(DEGs), ], scale = "row", cluster_rows = T,cluster_cols = T,color = colorRampPalette(c("black", "white", "firebrick3"))(50),
              show_rownames = F, show_colnames = T, cutree_rows=2)


gene_order = ph$tree_row$labels[ph$tree_row$order]
# ATAC associate with RNA-seq DEGs
DEGs_atac = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/DEGs/DEGs.reads.bed", stringsAsFactors = F)
DEGs_atac = DEGs_atac[match(gene_order, DEGs_atac$V4), ]
pheatmap(DEGs_atac[,c(7,8,5,9,6,10)], scale = "row", cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("blue", "black", "yellow"))(50),
         show_rownames = F, show_colnames = T, cutree_rows=2)

# normalize ATAC signal by median
boxplot(log2(DEGs_atac[,c(7,8,5,9,6,10)]+1))
summary(DEGs_atac[,c(7,8,5,9,6,10)])

DEGs_atac_norm = DEGs_atac
DEGs_atac_norm$V7 = DEGs_atac_norm$V7/66.0 
DEGs_atac_norm$V5 = DEGs_atac_norm$V5/46.0
DEGs_atac_norm$V6 = DEGs_atac_norm$V6/36.00 
DEGs_atac_norm$V8 = DEGs_atac_norm$V8/321.0
DEGs_atac_norm$V9 = DEGs_atac_norm$V9/144
DEGs_atac_norm$V10 = DEGs_atac_norm$V10/224 

pdf("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/figures/RNA.DEGs_ATAC_gene2K.heatmap.pdf", width = 5, height =7)
pheatmap(DEGs_atac_norm[,c(8,7,5,9,6,10)], scale = "row", cluster_rows = F,cluster_cols = F,color = colorRampPalette(c("purple", "black", "yellow"))(50),
         show_rownames = F, show_colnames = F, gaps_row=201, fontsize=20)

dev.off()


# ZNF gene cluster
ZNF = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA-seq_Replicate/ZNF.gene.cpm.txt", row.names = 1)
colnames(ZNF) = c("midlife_rep1", "old_rep1", "young_rep1", "old_rep2", "midlife_rep2", "young_rep2")
ZNF = ZNF[,c("young_rep1", "young_rep2", "midlife_rep1", "midlife_rep2", "old_rep1", "old_rep2")]
pheatmap(ZNF[rowMeans(ZNF) >= 10, ], scale = "row", cluster_cols = F)

boxplot(log2(ZNF+1))
refGene = read.table("/lustre/user/liclab/publicData/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/refGene.txt", sep="\t")
refGene = refGene[refGene$V13 %in% rownames(ZNF), c(3, 13)]

ZNF_chr19 = ZNF[unique(refGene[refGene$V3 == "chr19", 2]), ]
pheatmap(ZNF_chr19[rowMeans(ZNF_chr19) >= 3, ], scale = "row", cluster_cols = F)

boxplot(log2(ZNF_chr19[rowMeans(ZNF_chr19) >= 1, ] +1))




rep1 = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA/pipe_out/cufflinks/cuffnorm/genes.fpkm_table",  header = T, row.names = 1)
rep2 = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA_rep2/pipe_out/cufflinks/cuffnorm/genes.fpkm_table", header = T, row.names = 1)
znf_z2h2 = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA-seq_Replicate/ZNF-C2H2.txt", sep="\t")
rep1 = rep1[, -3]
expr = merge(rep1, rep2, by = "row.names")
rownames(expr) = expr$Row.names
expr = expr[,-1]
colnames(expr) = c("midlife_rep1", "old_rep1", "young_rep1", "old_rep2", "midlife_rep2", "young_rep2")
expr = expr[,c("young_rep1", "young_rep2", "midlife_rep1", "midlife_rep2", "old_rep1", "old_rep2")]

refGene = read.table("/lustre/user/liclab/publicData/igenomes/Homo_sapiens/UCSC/hg19/Annotation/Genes/refGene.txt", sep="\t")
refGene = refGene[refGene$V13 %in% znf_z2h2$V1, c(3, 13)]

expr_ZNF = na.omit(expr[refGene[refGene$V3 == "chr19", 2], ])
pheatmap(expr_ZNF[rowMeans(expr_ZNF) >= 1, c(2,4,6, 1,3,5)], scale = "row", cluster_cols = F)

pheatmap(expr_ZNF[expr_ZNF$old_rep1 < expr_ZNF$young_rep1 & expr_ZNF$old_rep2 < expr_ZNF$young_rep2, ], scale = "row", cluster_cols = F)

boxplot(log2(expr_ZNF$old_rep1/expr_ZNF$young_rep1), log2(expr_ZNF$old_rep2/expr_ZNF$young_rep2))

