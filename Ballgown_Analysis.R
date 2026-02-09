#Loading and preparing data in R for Ballgown
install.packages("ballgown")
library(ballgown)
library(dplyr)

#Load the Ballgown data(raw data)
bg <- ballgown(dataDir ="/home/asus/Desktop/embs/samples/sra/ballgown/", samplePattern ="SRR")
bg              

sample_ids <- c("SRR24491385", "SRR24491386", "SRR24491387", "SRR24491388",
                "SRR24491392", "SRR24491394", "SRR24491395", "SRR24491396",
                "SRR24491419", "SRR24491424", "SRR24491397", "SRR24491399",
                "SRR24491401", "SRR24491402", "SRR24491405", "SRR24491406",
                "SRR24491411", "SRR24491434", "SRR24491433", "SRR24491432")
group_labels <- c(rep("control", 10), rep("treatment", 10))
pData <- data.frame(sample=sample_ids, group=group_labels)

#View the sample data (for verifying control vs treatment labeling)
pData

texpr(bg) [1:3, 1:7]

install.packages("matrixStats")
library(matrixStats)

rowVars(texpr(bg)) >1

#filtered data
bg_filt= subset(bg, "rowVars(texpr(bg)) >1", genomesubset= TRUE)
bg_filt

pData(bg_filt) = pData 

dim(texpr(bg_filt))
dim(texpr(bg))

length(sampleNames(bg_filt))
list.dirs("/home/asus/Desktop/embs/samples/sra/ballgown", full.names = FALSE, recursive = FALSE)

#transcript level analysis
result_transcript= stattest(bg_filt, feature = "transcript", covariate = "group", getFC= TRUE, meas ="FPKM" )
result_transcript

head(result_transcript [order(result_transcript$pval), ])

#gene level analysis
result_gene= stattest(bg_filt, feature = "gene", covariate = "group", getFC= TRUE, meas= "FPKM")
result_gene

head(result_gene [order(result_gene$pval), ])

#add gene names and gene ids to transcript level
result_transcript= data.frame(geneNames= ballgown::geneNames(bg_filt), geneIDs= ballgown::geneIDs(bg_filt),result_transcript)
head(result_transcript)

dim(result_transcript)

dim(result_gene)

#filter the significant genes using p value
result_transcript=arrange(result_transcript, pval)
result_gene=arrange(result_gene, pval)
head(result_gene)

#write the results to csv file
write.csv(result_transcript, "transcripts_results.csv")
write.csv(result_gene, "genes_results.csv")



#additional part
#to add gene names and ids to gene level
#Read both CSV files
transcript_df <- read.csv("transcripts_results.csv")
gene_df <- read.csv("genes_results.csv")

#Create gene ID → gene name mapping from transcript-level data
gene_lookup <- unique(transcript_df[, c("geneIDs", "geneNames")])
colnames(gene_lookup) <- c("gene_id", "gene_name")

#Rename the gene ID column in gene_df to match for merging
colnames(gene_df)[colnames(gene_df) == "id"] <- "gene_id"

#Merge gene names into gene_df
gene_annotated <- merge(gene_df, gene_lookup, by = "gene_id", all.x = TRUE)

#Save result
write.csv(gene_annotated, "gene_level_with_gene_names.csv", row.names = FALSE)

#Check how many gene names were matched
sum(is.na(gene_annotated$gene_name))
head(gene_annotated)

#Prioritize non-missing gene names(to find duplicates)
gene_annotated_clean <- gene_annotated[gene_annotated$gene_name != ".", ]

#Remove duplicate gene_id entries, keeping first
gene_annotated_clean <- gene_annotated_clean[!duplicated(gene_annotated_clean$gene_id), ]

#Save cleaned version
write.csv(gene_annotated_clean, "gene_level_with_cleaned_names.csv", row.names = FALSE)

sum(is.na(gene_annotated_clean$gene_name))
head(gene_annotated_clean)

colnames(transcript_df)
colnames(gene_df)
colnames(gene_annotated_clean)

#reorder columns so that gene_name is second column
gene_annotated_clean <- gene_annotated_clean[, c("gene_id", "gene_name", setdiff(names(gene_annotated_clean), c("geneIDs", "geneNames")))]
head(gene_annotated_clean)

#Drop duplicated columns
gene_annotated_clean <- gene_annotated_clean[, !(colnames(gene_annotated_clean) %in% c("gene_id.1", "gene_name.1", "X"))]

#Reorder: gene_id, gene_name, then rest
gene_annotated_clean <- gene_annotated_clean[, c("gene_id", "gene_name", setdiff(colnames(gene_annotated_clean), c("gene_id", "gene_name")))]

#View result
head(gene_annotated_clean)



#final steps
#filter the significant transcripts by p value
final_transcript_list = subset(result_transcript, result_transcript$pval<=0.05)
head(final_transcript_list)

#filter the significant genes by p value
final_gene_list = subset(gene_annotated_clean, gene_annotated_clean$pval<=0.05)

#write the results to csv file
write.csv(final_transcript_list, "final_transcript_list.csv")
write.csv(final_gene_list, "final_gene_list.csv")

sum(result_gene$pval <= 0.05)
sum(result_transcript$pval <= 0.05)

colnames(gene_annotated_clean)
colnames(gene_lookup) <- c("gene_id", "gene_name")
colnames(gene_df)[colnames(gene_df) == "id"] <- "gene_id"


#VISUALIZATION
#BOX PLOT of raw data(bg)
fpkm = texpr(bg, meas = "FPKM")
fpkm
boxplot(fpkm)     #it won't show good, because it's a raw data. so we have to normalize it

#normalizing
boxplot(log2(fpkm+1), las= 2, xlab = "samples",ylab = "fpkm values", main = "log2(FPKM+1)Distribution")

#BOX PLOT of filtered data(bg_filt)
fpkm = texpr(bg_filt, meas = "FPKM")
boxplot(fpkm) 

#*************************
# BOX PLOT of raw data (bg)
  fpkm <- texpr(bg, meas = "FPKM")

# Adjust margins before plotting
par(mar = c(10, 4, 4, 2))   # Bottom, Left, Top, Right margins

# Normalizing and plotting with rotated labels
boxplot(
  log2(fpkm + 1),
  las = 2,                   # Vertical labels
  cex.axis = 0.8,            # Slightly smaller labels
  xlab = "Samples",
  ylab = "FPKM values",
  main = "log2(FPKM+1) Distribution"
)

# BOX PLOT of filtered data (bg_filt)
fpkm_filt <- texpr(bg_filt, meas = "FPKM")
par(mar = c(10, 4, 4, 2))
boxplot(
  log2(fpkm_filt + 1),
  las = 2,
  cex.axis = 0.8,
  xlab = "Samples",
  ylab = "FPKM values",
  main = "Filtered log2(FPKM+1) Distribution"
)


#normalizing (the filtered boxplot is good comparing to first one)
boxplot(log2(fpkm+1), las= 2, xlab = "samples",ylab = "fpkm values", main = "log2(FPKM+1)Distribution")

#************************
# Set consistent margin for long sample names
par(mar = c(10, 4, 4, 2))

# 1 Normalized RAW data
fpkm_raw <- texpr(bg, meas = "FPKM")
boxplot(
  log2(fpkm_raw + 1),
  las = 2,
  cex.axis = 0.8,
  xlab = "",
  ylab = "FPKM values",
  main = "Normalized Raw Data: log2(FPKM+1)"
)

# 2 Normalized FILTERED data
fpkm_filt <- texpr(bg_filt, meas = "FPKM")
par(mar = c(10, 4, 4, 2)) # ensure margins reset
boxplot(
  log2(fpkm_filt + 1),
  las = 2,
  cex.axis = 0.8,
  xlab = "",
  ylab = "FPKM values",
  main = "Normalized Filtered Data: log2(FPKM+1)"
)





ballgown::geneNames(bg_filt)[10]
plot(fpkm[10,]~as.factor(pData$group))

ballgown::geneNames(bg_filt)[20]
plot(fpkm[20,]~as.factor(pData$group))

#to get the gene name in 20th position
geneNames(bg_filt)[20]
results <- stattest(bg_filt, feature="gene", covariate="group", getFC=TRUE, meas="FPKM")
results[20, ]            #To get stats for this specific gene



#VOLCANO PLOT FOR DEGs
# Load libraries
library(ggplot2)
install.packages("readr")
library(readr)

# Load gene data
genes <- read_csv("final_gene_list.csv")

genes$log2fc <- log2(genes$fc)
genes$significance <- "Not Significant"

# Define significance criteria
genes$significance[genes$pval < 0.05 & genes$log2fc > 1] <- "Upregulated"
genes$significance[genes$pval < 0.05 & genes$log2fc < -1] <- "Downregulated"

# Plot
ggplot(genes, aes(x = log2fc, y = -log10(pval), color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs",
       x = "log2(Fold Change)",
       y = "-log10(p-value)") +
  theme(legend.title = element_blank())


install.packages("ggrepel")
library(ggrepel)

top_genes <- genes %>%
  dplyr::filter(pval < 0.05) %>%
  dplyr::arrange(pval) %>%
  dplyr::slice(1:10)

ggplot(genes, aes(x = log2fc, y = -log10(pval), color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_text_repel(data = top_genes, aes(label = gene_name), size = 3) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs",
       x = "log2(Fold Change)",
       y = "-log10(p-value)") +
  theme(legend.title = element_blank())



#VOLCANO PLOT FOR DETs
# Load libraries
library(ggplot2)
install.packages("readr")
library(readr)

# Load transcript data
transcripts <- read_csv("final_transcript_list.csv")

transcripts$log2fc <- log2(transcripts$fc)
transcripts$significance <- "Not Significant"

# Define significance criteria
transcripts$significance[transcripts$pval < 0.05 & transcripts$log2fc > 1] <- "Upregulated"
transcripts$significance[transcripts$pval < 0.05 & transcripts$log2fc < -1] <- "Downregulated"

# Plot
ggplot(transcripts, aes(x = log2fc, y = -log10(pval), color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red", 
                                "Downregulated" = "blue", 
                                "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of DEGs",
       x = "log2(Fold Change)",
       y = "-log10(p-value)") +
  theme(legend.title = element_blank())


install.packages("ggrepel")
library(ggrepel)

install.packages("dplyr")
library(dplyr)

top_transcripts <- transcripts %>%
  dplyr::filter(pval < 0.05) %>%
  dplyr::arrange(pval) %>%
  dplyr::slice(1:10)

ggplot(transcripts, aes(x = log2fc, y = -log10(pval), color = significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_text_repel(data = top_transcripts, aes(label = geneNames), size = 3) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Volcano Plot of DETs",
       x = "log2(Fold Change)",
       y = "-log10(p-value)") +
  theme(legend.title = element_blank())



#HEAT MAP for DEGs
library(dplyr)
install.packages("pheatmap")
library(pheatmap)
install.packages("dplyr")  
library(dplyr) 

# Load DEG list
deg_genes <- read.csv("final_gene_list.csv")

fpkm_mat <- texpr(bg_filt, meas = "FPKM")  
rownames(fpkm_mat) <- ballgown::geneNames(bg_filt)  # Set rownames to gene names

rownames(fpkm_mat) <- ballgown::geneIDs(bg_filt)
head(deg_genes_sig$gene_id)

intersect_ids <- intersect(rownames(fpkm_mat), deg_genes_sig$gene_id)
length(intersect_ids)

deg_fpkm_mat <- fpkm_mat[rownames(fpkm_mat) %in% deg_genes_sig$gene_id, ]
dim(deg_fpkm_mat)

# Log2 transform
deg_fpkm_mat_log2 <- log2(deg_fpkm_mat + 1)

# Z-score scaling by gene
scaled_deg_fpkm <- t(scale(t(deg_fpkm_mat_log2)))

annotation_col <- data.frame(Group = factor(group_labels))
rownames(annotation_col) <- colnames(scaled_deg_fpkm)

library(pheatmap)

pheatmap(scaled_deg_fpkm,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         fontsize_col = 8,
         main = "Heatmap of DEGs (log2 FPKM +1, scaled)")


#HEAT MAP showing top 50 significant DEGs
#Show only the top 50 most significant DEGs
# Sort DEG genes by p-value
deg_genes_sig <- deg_genes_sig[order(deg_genes_sig$pval), ]

# Select top 50 DEGs
top_deg_genes <- head(deg_genes_sig, 50)

# Get matching FPKM values
top_deg_fpkm_mat <- fpkm_mat[top_deg_genes$gene_id, ]

# Assign gene names as row names
rownames(top_deg_fpkm_mat) <- top_deg_genes$gene_name

# Log2 transform
top_deg_fpkm_mat_log2 <- log2(top_deg_fpkm_mat + 1)

# Z-score scaling by gene
scaled_top_deg_fpkm <- t(scale(t(top_deg_fpkm_mat_log2)))

# Sample annotation
annotation_col <- data.frame(Group = factor(group_labels))
rownames(annotation_col) <- colnames(scaled_top_deg_fpkm)

# Plot the top 50 DEG heatmap
pheatmap(scaled_top_deg_fpkm,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 8,         # Slightly larger
         fontsize_col = 8,
         main = "Top 50 DEGs (log2 FPKM +1, scaled)")



#DAVID
#Extract Gene Symbols from final_gene_list.csv
# Load your significant DEGs
deg_genes <- read.csv("final_gene_list.csv")

# Extract the gene symbols
gene_symbols <- deg_genes$gene_name

# Check the first few gene names
head(gene_symbols)

# Write gene symbols to a plain text file (one per line)
gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != "."]
write.table(gene_symbols, file = "deg_gene_symbols.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)



#***********************************************


#Step-by-step to count up- and downregulated genes
# Load the final gene list
deg_genes <- read.csv("final_gene_list.csv")

# Ensure log2fc and significance columns exist
deg_genes$log2fc <- log2(deg_genes$fc)
deg_genes$significance <- "Not Significant"
deg_genes$significance[deg_genes$pval < 0.05 & deg_genes$log2fc > 1] <- "Upregulated"
deg_genes$significance[deg_genes$pval < 0.05 & deg_genes$log2fc < -1] <- "Downregulated"

# Count the number of up- and downregulated genes
table(deg_genes$significance)

upregulated_count <- sum(deg_genes$significance == "Upregulated")
downregulated_count <- sum(deg_genes$significance == "Downregulated")
not_significant_count <- sum(deg_genes$significance == "Not Significant")

cat("Upregulated:", upregulated_count, "\n")
cat("Downregulated:", downregulated_count, "\n")
cat("Not Significant:", not_significant_count, "\n")

# Subset upregulated and downregulated genes
library(dplyr)

# Read in your final gene list
genes <- read.csv("final_gene_list.csv")

# Check the first few rows
head(genes)

# Calculate log2 fold change
genes$log2fc <- log2(genes$fc)

#Add log2fc and significance
# Add significance category
genes$significance <- "Not Significant"
genes$significance[genes$pval < 0.05 & genes$log2fc > 1] <- "Upregulated"
genes$significance[genes$pval < 0.05 & genes$log2fc < -1] <- "Downregulated"

#Check if 'significance' column was created
table(genes$significance)
colnames(genes)

#Group by and count
gene_counts <- genes %>%
  group_by(significance) %>%
  summarise(count = n())

#Bar plot
ggplot(gene_counts, aes(x = significance, y = count, fill = significance)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Differentially Expressed Genes",
       x = "Category",
       y = "Number of Genes") +
  theme(legend.position = "none")


#for transcript
transcripts <- read.csv("final_transcript_list.csv")

transcripts$log2fc <- log2(transcripts$fc)
transcripts$significance <- "Not Significant"
transcripts$significance[transcripts$pval < 0.05 & transcripts$log2fc > 1] <- "Upregulated"
transcripts$significance[transcripts$pval < 0.05 & transcripts$log2fc < -1] <- "Downregulated"

transcript_counts <- transcripts %>%
  group_by(significance) %>%
  summarise(count = n())

ggplot(transcript_counts, aes(x = significance, y = count, fill = significance)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "grey")) +
  theme_minimal() +
  labs(title = "Differentially Expressed Transcripts",
       x = "Category",
       y = "Number of Transcripts") +
  theme(legend.position = "none")

#Create DEG Summary Table
library(dplyr)

genes <- read.csv("final_gene_list.csv")
genes$log2fc <- log2(genes$fc)
genes$padj <- p.adjust(genes$pval, method = "BH")

#Annotate DEG Category
genes$category <- "Not Significant"
genes$category[genes$padj < 0.05 & genes$log2fc > 1] <- "Upregulated"
genes$category[genes$padj < 0.05 & genes$log2fc < -1] <- "Downregulated"

#Filter for DEGs Only
degs <- genes %>%
  filter(category != "Not Significant")

#Select Relevant Columns
deg_summary <- degs %>%
  select(gene_name, log2fc, pval, padj, category)

colnames(genes)

write.csv(deg_summary, "DEG_summary_table.csv", row.names = FALSE)

