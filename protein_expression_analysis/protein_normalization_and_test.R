library("DEP")
library("dplyr")
library("readxl")
# Read the proteomics dataset from the specified Excel sheet
Dataset <- read_excel("./input/proteome_for_DEP.xlsx", sheet='Liv_WT_vs_ob')
data_unique <- make_unique(Dataset, "Gene.Name", "Protein.ID", delim = ";")
iBAQ_columns <- grep("_", colnames(data_unique))
data_se_parsed <- make_se_parse(data_unique, iBAQ_columns)
data_filt2 <- filter_missval(data_se_parsed, thr = 2) # Filter out proteins with too many missing values (keep proteins with >= thr valid values per group)
data_norm <- normalize_vsn(data_filt2) # Normalize intensity values using variance-stabilizing normalization (VSN)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01) # Impute missing values using the MinProb method (shifted distribution), q = 0.01
data_diff <- test_diff(data_imp, type = "control", control = "WT_") # Perform differential expression testing between groups
dep <- add_rejections(data_diff, alpha = 0.05)
# Extract the table of statistical test results
data_results <- get_results(dep)
write.table(data_results, "./output/WT_ob_proteome_comparison.txt", quote=F)
write.table(data_norm@assays@data@listData, "./output/proteome_after_norm_WT_ob.txt", quote=F)
write.table(data_norm@elementMetadata@listData, "./output/proteome_after_norm_protein_name_WT_ob.txt", quote=F)
