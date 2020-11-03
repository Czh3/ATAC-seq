# super enhancer
SE = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/SE.reads_reform.bed", row.names = 1)
colnames(SE) = c( "midlife.rep2", "old.rep2", "young.rep2", "young.rep1", "midlife.rep1", "old.rep1")
pheatmap(SE[,1:3], scale = "row", show_rownames = F)
pheatmap(SE[,4:6], scale = "row", show_rownames = F)
SE$FC_rep1 = log2(SE$old.rep1/SE$young.rep1)
SE$FC_rep2 = log2(SE$old.rep2/SE$young.rep1)
boxplot(SE$FC_rep1 , SE$FC_rep2)


se_gene = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/ENCODE/SE.gene.bed")
se_gene = unique(se_gene$V6)

gene_expr = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/RNA-seq_Replicate/RNA.cpm.txt", header = T, row.names = 1)

se_gene_expr = na.omit(gene_expr[se_gene, ])
boxplot(log2(se_gene_expr+1))
pheatmap(se_gene_expr, scale = "row", show_rownames = F)


# RAS
SE_RAS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/OIS/ATAC_new/peaks/merge/OIS.SE.reads.bed")
boxplot(log10(SE_RAS[,6:11]))
pheatmap(SE_RAS[6:11], scale = "row", show_rownames = F)
