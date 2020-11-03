# repeat masker
setwd("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/Rscript")

Constant = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/repeatMask/contantPeaks.repeatMask.bed")
IARs = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/repeatMask/edgeR.increase.ATAC.peaks.repeatMask.bed")
DARs = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/repeatMask/edgeR.decrease.ATAC.peaks.repeatMask.bed")

rp = data.frame(rbind(cbind(Constant$V10, "Constant"),
           cbind(IARs$V15, "IARs"),
           cbind(DARs$V15, "DARs")))

colnames(rp) = c("repeats", "samples")

ggplot(rp, aes(factor(repeats), fill = samples)) + 
  geom_bar(position=position_dodge()) +
  theme_classic(base_size = 25) + xlab("") + ylab("percent") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.line = element_line(size = 1, colour = "black"))


rp1 = data.frame(rbind(cbind(Constant$V11, "Constant"),
                      cbind(IARs$V16, "IARs"),
                      cbind(DARs$V16, "DARs")))

colnames(rp1) = c("repeats", "samples")

ggplot(rp1, aes(factor(repeats), fill = samples)) + 
  geom_bar(position=position_dodge()) +
  theme_classic(base_size = 25) + xlab("") + ylab("percent") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.line = element_line(size = 1, colour = "black"))
