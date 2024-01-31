pacman::p_load(DESeq2, readxl, dplyr, ReportingTools, ggplot2, EnhancedVolcano,
               writexl, tibble, rtracklayer, rlang, flextable)

format_scientific <- function(x) {
  formatC(x, format = "e", digits = 2)
}

leaf_counts <- read_excel("PATH_HERE/transcripts.readcount.annot.xlsx") |> 
  rename("Leaf1"=Input_1,
         "Leaf2"=Input_2,
         "Leaf3"=Input_3) |> 
  select(-gene_name, -gene_description)
phloem_counts <- read_excel("PATH_HERE/transcripts.readcount.annot.xlsx")|> 
  select(-gene_name, -gene_description)


# Getting lncRNAS gene ids
leafCP <- read_excel("PATH_HERE/LncRNA_Filter/Coding_potential_filter.result.xlsx")
phloemCP <- read_excel("PATH_HERE/Coding_potential_filter.result.xlsx")

leaf_lncRNAs <- leafCP |> filter(CNCI=='noncoding' & CPC=='noncoding' & PFAM=='noncoding')
leaf_lncRNAs_TCONS <- leaf_lncRNAs |> 
  select(transcript_id) |> 
  unique() #removing duplicates

phloem_lncRNAs <- phloemCP |> filter(CNCI=='noncoding' & CPC=='noncoding' & PFAM=='noncoding')
phloem_lncRNAs_TCONS <- phloem_lncRNAs |> 
  select(transcript_id) |> 
  unique()


# Getting lncRNAs counts and filtering out counts<5
leaf_tracking_lncRNA_IDs <- merge(leaf_counts, leaf_lncRNAs_TCONS) |> 
  filter(!Leaf1<5 & !Leaf2<5 & !Leaf3<5) |> 
  relocate(gene_id) |> 
  rename(Leaf_id = transcript_id)

leaf_count_sums <- aggregate(. ~ gene_id, data = leaf_tracking_lncRNA_IDs [-2], FUN = sum)

phloem_tracking_lncRNA_IDs <- merge(phloem_counts, phloem_lncRNAs_TCONS) |> 
  filter(!Phloem1<5 & !Phloem2<5 & !Phloem3<5) |> 
  relocate(gene_id) |> 
  rename(Phloem_id = transcript_id)


# Merging and summing isoforms counts
phloem_count_sums <- aggregate(. ~ gene_id, data = phloem_tracking_lncRNA_IDs[-2], FUN = sum)


# Running DESeq2
samples_counts <-  merge(phloem_count_sums, leaf_count_sums)|> 
  filter(!grepl(paste0("^", "XLOC_"), gene_id)) #removing XLOC ids (they come from cufflinks, most likely weren't mapped to an annotated gene)
row.names(samples_counts) <-  samples_counts$gene_id
samples_counts <- samples_counts[-1]

samples_info <- tibble(Sample=c("Phloem1", "Phloem2", "Phloem3", "Leaf1", "Leaf2", "Leaf3"),
                       Tissue=c("phloem", "phloem", "phloem", "leaf", "leaf", "leaf"))
samples_info$Tissue <- factor(samples_info$Tissue)

sample_dds <- DESeqDataSetFromMatrix(countData = samples_counts,
                                     colData = samples_info,
                                     design = ~Tissue)
sample_dds <- DESeq(sample_dds)

sample_results <- results(sample_dds,contrast = c("Tissue", "phloem", "leaf"))


# Saving results in a table
results_table <- sample_results |> data.frame() |>
  rownames_to_column(var = "gene_id") |> 
  arrange(padj) |> 
  mutate_at(vars(padj), ~(signif(., 3))) |> 
  mutate_at(vars(baseMean, log2FoldChange, lfcSE, stat), ~(round(., 2))) |> 
  select(-stat, -pvalue)

#
DE_counts <- merge(results_table, phloem_count_sums)
DE_counts <-  merge(DE_counts, leaf_count_sums)
write_xlsx(DE_counts, path = "DESeq2/DESeq2_PhloemvsLeaf_lncRNAs_wCounts.xlsx")
#

write_xlsx(results_table, path = "DESeq2/DESeq2_PhloemvsLeaf_lncRNAs.xlsx")


# Plotting
## Volcano
keyvals <- ifelse(
  sample_results$log2FoldChange < -1 & sample_results$padj < 1e-5, 'red',
  ifelse(sample_results$log2FoldChange > 1 & sample_results$padj < 1e-5, 'green',
         'grey'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'green'] <- 'high'
names(keyvals)[keyvals == 'grey'] <- 'mid'
names(keyvals)[keyvals == 'red'] <- 'low'

lncVolcano <- EnhancedVolcano(sample_results,
                              lab = rownames(sample_results),
                              x = 'log2FoldChange',
                              y = 'padj',
                              title = '',
                              subtitle = NULL,
                              axisLabSize = 18,
                              gridlines.major = F,
                              gridlines.minor = F,
                              colCustom = keyvals,
                              selectLab = NA,
                              xlab = noquote(expression('log'[2]*'FoldChange')),
                              ylab = noquote(expression('-log'[10]*'(padj)')),
                              labSize = 4,
                              pointSize = 4,
                              legendPosition = 'right',
                              legendLabSize = 18,
                              legendIconSize = 8,
                              drawConnectors = TRUE,
                              widthConnectors = .5,
                              colConnectors = 'black')+
  coord_fixed(ratio = .15)


ggsave(filename = "DESeq2/DESeq2_lncRNAs_PhloemvsLeaf_VolcanoPlot.png",
       plot = lncVolcano,
       width = 10, height = 8,
       dpi = 600)


# Getting significant results in a table
flextable(results_table[results_table$padj<=1e-5,]) |> 
  fontsize(size = 12, part = "header") |> 
  bold(i = NULL, j = NULL, bold = TRUE, part = "header") |> 
  align_text_col(align = "center", header = T, footer = T) |> 
  align_nottext_col(align = "center") |> 
  set_formatter(padj = format_scientific) |> 
  save_as_image(path = "DESeq2/DESeq2__PhloemvsLeaf_lncRNAs.png")


# Getting down-regulated
down <- results_table[results_table$log2FoldChange<0 & results_table$padj<=1e-5, ] #34
flextable(down) |> 
  fontsize(size = 12, part = "header") |> 
  bold(i = NULL, j = NULL, bold = TRUE, part = "header") |> 
  align_text_col(align = "center", header = T, footer = T) |> 
  align_nottext_col(align = "center") |> 
  set_formatter(padj = format_scientific) |> 
  save_as_image(path = "DESeq2/DESeq2_PhloemvsLeaf_DOWN_lncRNAs.png")


# Getting up-regulated
up <- results_table[results_table$log2FoldChange>0 & results_table$padj<=1e-5, ] #33
flextable(up) |> 
  fontsize(size = 12, part = "header") |> 
  bold(i = NULL, j = NULL, bold = TRUE, part = "header") |> 
  align_text_col(align = "center", header = T, footer = T) |> 
  align_nottext_col(align = "center") |> 
  set_formatter(padj = format_scientific) |> 
  save_as_image(path = "DESeq2/DESeq2_PhloemvsLeaf_UP_lncRNAs.png")


###################################################################################
# Getting the gene_names for GO term analysis
annot <- readGFF("PATH_HERE/novel_lncRNA.gtf") |> 
  select(gene_id, gene_name) |> 
  unique()

## For all significant differential expressed
TopDE_gene_names <- merge(results_table, annot) |> 
  filter(padj<1e-5)
write_xlsx(TopDE_gene_names, path = "PATH_HERE/DESeq2_PhloemvsLeaf_ToplncRNAs_gene_names.xlsx")

## For all significant down regulated
merge(down, annot) |> 
write_xlsx(path = "PATH_HERE/DESeq2_PhloemvsLeaf_DOWNlncRNAs_gene_names.xlsx")

## For all significant up regulated
merge(up, annot) |> 
  write_xlsx(path = "PATH_HERE/DESeq2_PhloemvsLeaf_UPlncRNAs_gene_names.xlsx")
