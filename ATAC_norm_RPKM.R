# atac-seq normalized by RPKM

atac = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/merged_peaks.reads_reform.bed",row.names = 1)


colnames(atac) = c( "midlife.rep2", "old.rep2", "young.rep2", "young.rep1", "midlife.rep1", "old.rep1")
condition = factor(c("midlife", "old", "young", "young", "midlife", "old"), levels = c("young", "midlife", "old"))


library(DESeq)
cds = newCountDataSet( atac, condition )
cds = estimateSizeFactors( cds )
sizeFactors( cds )
count.table.norm = counts( cds, normalized=TRUE )

# pca
condition.info = data.frame("condition" = condition)
row.names(condition.info) = c( "midlife.rep2", "old.rep2", "young.rep2", "young.rep1", "midlife.rep1", "old.rep1")
pca = prcomp(t(log10(count.table.norm+1)), scale = T)
pca = as.data.frame(pca$x)

pca = cbind(condition, pca)
ggplot(pca, aes(PC1, PC2)) +
  geom_point(aes(color = condition), size=4)+
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_line(size = 1, colour = "black"),
        axis.text = element_text(colour = "black")) +
  scale_colour_brewer(palette="Dark2") 


cds = estimateDispersions( cds )

d = dist(t(log10(count.table.norm+1)))
h = hclust(d)
plot(h)

res = nbinomTest( cds, "young", "old" )

resSig = res[ res$padj < 0.05, ]
resSig = na.omit(resSig)
resSig <- resSig[ order(resSig$pval), ]


y <- DGEList(counts=atac,group=condition)
y <- calcNormFactors(y)
# normlized counts
logcpm <- cpm(y, prior.count=1, log=TRUE)
peaks_cpm <- cpm(y, prior.count=1, log=F)

design <- model.matrix(~condition)
y <- estimateDisp(y,design)

fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=3)
sigRegions = topTags(lrt, n=nrow(logcpm))$table
DARs = sigRegions[sigRegions$FDR <= 0.05, ]
dim(DARs)
