rm (list = ls())# Remove all objects from the current R environment

sample_type = "WT_00min_vs_OB_00min"

pair <- c("WT_00min","OB_00min")

library("edgeR")
# Read the gene expression count matrix (TSV format), with the first column as row names (gene IDs)
count <- read.table(paste("./input/Liver_00_WT11OB12_GeneExpressionMatrix.tsv", sep = ""), sep = "\t", header = T, row.names = 1)
count <- as.matrix(count)
# Assign column names to the count matrix (sample names)
col_name_list <- c("WT_00min_1", "WT_00min_2", "WT_00min_3", "WT_00min_4", "WT_00min_5", "WT_00min_6","WT_00min_7", "WT_00min_8", "WT_00min_9","WT_00min_10", "WT_00min_11",
                   "OB_00min_1", "OB_00min_2", "OB_00min_3","OB_00min_4", "OB_00min_5", "OB_00min_6","OB_00min_7", "OB_00min_8", "OB_00min_9","OB_00min_10", "OB_00min_11", "OB_00min_12")
colnames(count) <- col_name_list
# Define group labels for each sample
group_all <- factor(c("WT_00min", "WT_00min", "WT_00min","WT_00min", "WT_00min", "WT_00min","WT_00min", "WT_00min", "WT_00min","WT_00min", "WT_00min",
                      "OB_00min", "OB_00min", "0B_00min","OB_00min", "OB_00min", "OB_00min","OB_00min", "OB_00min", "OB_00min","OB_00min", "OB_00min", "OB_00min"))
# Filter out lowly expressed genes 
keep.exprs <- filterByExpr(count, group = group_all, keep.lib.sizes = FALSE)
count <- count[keep.exprs, ]
dim(count)
d <- DGEList(counts = count, group = group_all) # Create a DGEList object for edgeR
d <- calcNormFactors(d) # Calculate normalization factors for library size differences
d <- estimateCommonDisp(d) # Estimate the common dispersion
d <- estimateTagwiseDisp(d) # Estimate gene-specific (tagwise) dispersion

result <- exactTest(d, pair = pair)# Perform an exact test for differential expression between the two groups
# Extract the ranked table of results for all genes
table <- as.data.frame(topTags(result, n = nrow(count))) 
write.table(table, file = paste("./output/", sample_type, "_DEG.txt", sep = ""), col.names = T, row.names = T, sep = "\t")
# Create an MA plot highlighting the significant DEGs
is.DEG <- as.logical(table$FDR < 0.1)
DEG.names <- rownames(table)[is.DEG] 
png(paste("./output/", sample_type, "_MA_plot_Gene_FDR0.1.png", sep = ""))
plotSmear(result, de.tags = DEG.names, cex=0.3)
dev.off()
