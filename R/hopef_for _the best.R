if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("BiocManager")
BiocManager::install(version = "3.22")


core_packages<- c("tidyverse","ggpubr","ggsci","cowplot","patchwork")
install.packages(core_packages)

bioconductor_packages <- c("tximport", "DESeq2", "EnsDb.Hsapiens.v86","EnhancedVolcano")
BiocManager::install(bioconductor_packages)

install.packages("RColorBrewer")

BiocManager::install("AnnotationHub")
library(AnnotationHub)

BiocManager::install("restfulr")
library(EnsDb.Hsapiens.v86)

library(GenomicFeatures)
library(EnsDb.Hsapiens.v86)





install.packages("BiocManager")
BiocManager::install(
  "EnsDb.Hsapiens.v86",
  ask = FALSE,
  update = FALSE
)



samples <- list.dirs('input/salmon_out', recursive = FALSE, full.names = FALSE)
samples

quant_files <- file.path('input/salmon_out', samples, 'quant.sf')
quant_files


names(quant_files) <- samples
print(quant_files)

file.exists(quant_files)


col_data <- data.frame(
  row.names = samples,
  condition = factor(
    c(
      "HMGB2_CRISPR", "HMGB2_CRISPR", "HMGB2_CRISPR",
      "NT_CRISPR", "NT_CRISPR", "NT_CRISPR"
    )
  )
)

col_data

col_data$condition <- factor(col_data$condition)


columns(EnsDb.Hsapiens.v86)
keys(EnsDb.Hsapiens.v86)


tx2gene <- AnnotationDbi::select(EnsDb.Hsapiens.v86, 
                                 keys = keys(EnsDb.Hsapiens.v86),
                                 columns = c('TXID', 'SYMBOL'))

# check tx2gene 
head(tx2gene)
tx2gene <- dplyr::select(tx2gene, -GENEID)
head(tx2gene)

# Compile the tximport counts object and make DESeq dataset
# Get tximport counts object
library(tximport)

txi <- tximport(files = quant_files, 
                type = 'salmon',
                tx2gene = tx2gene,
                ignoreTxVersion = TRUE)

View(txi)
class(txi)


raw_counts  <- txi$counts
write.csv(raw_counts, "input/GSE288133_raw.csv", row.names = FALSE)

norm_counts <- txi$abundance

write.csv(norm_counts, "input/GSE288133_tpm.csv", row.names = FALSE)
txi_counts <- as.data.frame(txi)
names(txi_counts)
library(DESeq2)

dds <- DESeqDataSetFromTximport(txi = txi,
                                colData = col_data,
                                design = ~condition)



dds_rlog <- rlog(dds,blind = FALSE)

graphics.off()
install.packages("ggplot2")

plotPCA(dds_rlog, ntop = 500) +
theme_bw()

vsd <- vst(dds)
plotPCA(vsd)



dds <- DESeq(dds)


resdf <- results(dds)
View(resdf)

plotMA(resdf)

resdf <- as.data.frame(resdf)
resdf$gene <- rownames(resdf)
rownames(resdf) <- NULL

write.csv(resdf,file = "output/DESeq2_resilts.csv",row.names = FALSE)
saveRDS(resdf,file = "output/DESeq2_results.rds")





 