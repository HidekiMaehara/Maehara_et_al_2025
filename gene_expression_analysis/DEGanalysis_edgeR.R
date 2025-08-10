rm (list = ls())

sample_type = "WT_00min_vs_OB_00min"

pair <- c("WT_00min","OB_00min")

library("edgeR")
count <- read.table(paste("./input/Liver_00_WT11OB12_GeneExpressionMatrix.tsv", sep = ""), sep = "\t", header = T, row.names = 1)
count <- as.matrix(count)

col_name_list <- c("WT_00min_1", "WT_00min_2", "WT_00min_3", "WT_00min_4", "WT_00min_5", "WT_00min_6","WT_00min_7", "WT_00min_8", "WT_00min_9","WT_00min_10", "WT_00min_11",
                   "OB_00min_1", "OB_00min_2", "OB_00min_3","OB_00min_4", "OB_00min_5", "OB_00min_6","OB_00min_7", "OB_00min_8", "OB_00min_9","OB_00min_10", "OB_00min_11", "OB_00min_12")
colnames(count) <- col_name_list

group_all <- factor(c("WT_00min", "WT_00min", "WT_00min","WT_00min", "WT_00min", "WT_00min","WT_00min", "WT_00min", "WT_00min","WT_00min", "WT_00min",
                      "OB_00min", "OB_00min", "0B_00min","OB_00min", "OB_00min", "OB_00min","OB_00min", "OB_00min", "OB_00min","OB_00min", "OB_00min", "OB_00min"))

keep.exprs <- filterByExpr(count, group = group_all, keep.lib.sizes = FALSE)
count <- count[keep.exprs, ]
dim(count)
d <- DGEList(counts = count, group = group_all)
d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)

result <- exactTest(d, pair = pair)

table <- as.data.frame(topTags(result, n = nrow(count))) #toptagsで全行指定して???ータフレー???に変換
write.table(table, file = paste("./output/", sample_type, "_DEG.txt", sep = ""), col.names = T, row.names = T, sep = "\t")
is.DEG <- as.logical(table$FDR < 0.1) #logical クラス(true、false)に判定付きで変換???
DEG.names <- rownames(table)[is.DEG] #true???けDEG.namesに収???
png(paste("./output/", sample_type, "_MA_plot_Gene_FDR0.1.png", sep = ""))
plotSmear(result, de.tags = DEG.names, cex=0.3)
dev.off()
