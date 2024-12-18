library(gRodon)
library(Biostrings)

# seeting seed!
set.seed(1)

args <- commandArgs(trailingOnly = TRUE)
sample <- args[1]
input <- args[2]
output <- args[3]

genes <- readDNAStringSet(input)
highly_expressed <- grepl("ribosomal (subunit )?protein", names(genes), ignore.case = T)

res_10 <- predictGrowth(genes, highly_expressed, temperature = 10)
df_res_10 <- as.data.frame(res_10)
rownames(df_res_10) <- c(sample)

res <- predictGrowth(genes, highly_expressed, n_le = 1000)
df_res <- as.data.frame(res)
rownames(df_res) <- c(sample)

write.csv(df_res_10, output, row.names=TRUE)
