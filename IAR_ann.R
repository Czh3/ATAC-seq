# anno IAR/DAR

IARs_CS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/edgeR.increase.ATAC.peaks.Chromstate.bed")
IARs_CS = as.data.frame(summary(IARs_CS$V12)) / nrow(IARs_CS)
DARs_CS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/edgeR.decrease.ATAC.peaks.Chromstate.bed")
DARs_CS = as.data.frame(summary(DARs_CS$V12)) / nrow(DARs_CS)
AllPeaks_CS = read.table("/lustre/user/liclab/zhangc/Taolab/xiaohuangli/ATAC/rep2/peaks_DE/merged_peaks.shuffle.Chromstate.bed")
AllPeaks_CS = as.data.frame(summary(AllPeaks_CS$V7)) / nrow(AllPeaks_CS)


chromState = cbind(IARs_CS, DARs_CS, AllPeaks_CS)
chromState$IAR = log2(chromState$`summary(IARs_CS$V12)` / chromState$`summary(AllPeaks_CS$V7)`)
chromState$DAR = log2(chromState$`summary(DARs_CS$V12)` / chromState$`summary(AllPeaks_CS$V7)`)

chromState = chromState[,c(4,5)]
chromState$chromState = rownames(chromState)

chromState.m = melt(chromState)
chromState.m$chromState = factor(chromState.m$chromState,
                                 levels = rev(c("15_Repetitive/CNV","14_Repetitive/CNV", "13_Heterochrom/lo", "12_Repressed",
                                                "11_Weak_Txn", "10_Txn_Elongation", "9_Txn_Transition", "8_Insulator", "7_Weak_Enhancer",
                                                "6_Weak_Enhancer", "5_Strong_Enhancer", "4_Strong_Enhancer", "3_Poised_Promoter", "2_Weak_Promoter", "1_Active_Promoter")))



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
ggsave("../figures/IAR_DAR_chromStat.pdf", device = "pdf", width = 8, height = 5)

#
chromState = cbind(IARs_CS, DARs_CS, AllPeaks_CS)
colnames(chromState) = c("IAR","DAR","random")
chromState$chromState = rownames(chromState)

chromState.m = melt(chromState)

chromState.m$chromState = factor(chromState.m$chromState,
                                 levels = rev(c("15_Repetitive/CNV","14_Repetitive/CNV", "13_Heterochrom/lo", "12_Repressed",
                                                "11_Weak_Txn", "10_Txn_Elongation", "9_Txn_Transition", "8_Insulator", "7_Weak_Enhancer",
                                                "6_Weak_Enhancer", "5_Strong_Enhancer", "4_Strong_Enhancer", "3_Poised_Promoter", "2_Weak_Promoter", "1_Active_Promoter")))


ggplot(data=chromState.m, aes(x=chromState , y=value, fill=variable)) +
  geom_bar(stat="identity", position=position_dodge()) +
  #coord_flip() +
  #geom_hline(yintercept=0, size = 1, color = "black") +
  theme_classic(base_size = 20) + ylab("Percent of peaks") + xlab("Chrom State") +
  scale_color_brewer(palette="Set1") +
  scale_fill_brewer(palette="Set1") +
  theme(axis.line = element_line(size = 1, colour = "black"),
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        axis.text = element_text(colour = "black"),
        #axis.line.y = element_line(size = 0, colour = "black"),
        legend.position = "top")
ggsave("../figures/IAR_DAR_chromStat_percent.pdf", device = "pdf", width = 10, height = 5)


