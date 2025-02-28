if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(limma)

library(affy)

getwd()

# working directory 
setwd("D:/DEG analysis_oisharja/miRNA profile")

# read cell files
files <- list.files(pattern = "\\.CEL$", ignore.case = TRUE)
rawData <- ReadAffy(filenames = files, celfile.path = "D:/DEG analysis_oisharja/miRNA profile")

# RMA bg correction and normalization
eset <- rma(rawData)

# phenotype data of expressionset matrix 
head(pData(eset)) # manually assigning groups 

# grouping

group <- rep("Tumor", 25) #default grouping of all samples 

group[13:18] <- "Normal" # assigning 13-18th sample as normal 

group <- factor(group, levels = c("Tumor", "Normal")) # converting vector to factor

pData(eset)$group <- group # adding group to phenotype data of expressionset

head(pData(eset), 25) # verification 

# design matrix 
design <- model.matrix(~0 + group)

colnames(design) <- levels(group) # assigning colnames according to levels

# DEG 

# fit the linear model to the expression matrix 
fit <- lmFit(exprs(eset), design) 

# defining contrast: Tumor vs Normal 
contraxt.matrix <- makeContrasts(Tumor - Normal, levels = design)
fit2 <- contrasts.fit(fit, contraxt.matrix)

# empirical bayes moderation 
fit2 <- eBayes(fit2)

# extract DE miRNA 
results <- topTable(fit2, adjust.method = "BH", number = Inf)
head(results)

# significant genes
sig_results <- results[results$adj.P.Val < 0.05 & abs(results$logFC) > 1,]

head(sig_results)

# export 
write.csv(sig_results, file = "significant_miRNA.csv", row.names = TRUE)

# visualization 

# volcano plot 
library(ggplot2)

# convert results as dataframe and make a column with miRNA names
df_results <- as.data.frame(results)
df_results$miRNA <- rownames(df_results)

# marking significant miRNA 
df_results$Significant <- ifelse(df_results$adj.P.Val < 0.05 & abs(df_results$logFC) > 1, "Significant", "Not Significant")

# plot
ggplot(df_results, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(aes(color = Significant), size = 1) +
  scale_color_manual(values = c("Not significant" = "grey", "Significant" = "red")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Differentially Expressed miRNA",
    x = "Log2 Fold Change", 
    y = "-Log10(Adj.P.Val)"
  )


# heatmap

library(pheatmap)

# identifying significxant miRNAs
sig_miRNA <- rownames(results)[results$adj.P.Val < 0.05 & abs(results$logFC) > 1]

# subset expression matrix
sig_expr <- exprs(eset)[sig_miRNA, ]

# heatmap
pheatmap(sig_expr,
         scale = "row",
         annotation_col = pData(eset)[, "group", drop = FALSE],
         main = "Heatmap of Significant miRNAs")

# heatmap with top genes 

# top 50 by adjusted p-value
top_sig_miRNA <- head(rownames(results[order(results$adj.P.Val), ]), 50)

sig_expr_50 <- exprs(eset)[top_sig_miRNA, ]

pheatmap(sig_expr_50,
         scale = "row",
         annotation_col = pData(eset)[, "group", drop = FALSE],
         main = "Heatmap of Top 50 miRNAs")


